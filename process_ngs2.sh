#!/bin/bash
#SBATCH --job-name=process_ngs          # Job name
#SBATCH --output=condor.log             # Log for standard output
#SBATCH --error=condor.err              # Log for standard error
#SBATCH --cpus-per-task=40              # Adjust based on your ${cpus} variable
#SBATCH --mem=900G                       # Adjust based on your ${memory} variable
#SBATCH --partition=highcpu
# Set TMPDIR to /scratch to store temporary files there
export TMPDIR=/scratch/$USER/$SLURM_JOB_ID
#SBATCH --tmp=400G
#SBATCH --time=30:00:00                 # Maximum runtime (adjust as needed)
# Set TMPDIR to /scratch to store temporary files there
set -x
# Create the directory if it doesn't exist

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <params_file.csv>"
  exit 1
fi

params_file=$1

# Check if the params file exists
if [ ! -f "$params_file" ]; then
  echo "Error: Params file '$params_file' not found."
  exit 1
fi

# Get parameters from params file
path=$(awk -F, '/^data_filepath/ {print $2}' "$params_file" | tr -d '\r')
pear=$(awk -F, '/^pear_filepath/ {print $2}' "$params_file" | tr -d '\r')
pear_overlap=$(awk -F, '/^pear_overlap/ {print $2}' "$params_file" | tr -d '\r')
pear_stattest=$(awk -F, '/^pear_stattest/ {print $2}' "$params_file" | tr -d '\r')
pear_pvalue=$(awk -F, '/^pear_pvalue/ {print $2}' "$params_file" | tr -d '\r')
merge=$(awk -F, '/^merge/ {print $2}' "$params_file" | tr -d '\r')
compiler=$(awk -F, '/^compiler_filepath/ {print $2}' "$params_file" | tr -d '\r')
cpus=$(awk -F, '/^cpus/ {print $2}' "$params_file" | tr -d '\r')
memory=$(awk -F, '/^memory/ {print $2}' "$params_file" | tr -d '\r')
disk=$(awk -F, '/^disk/ {print $2}' "$params_file" | tr -d '\r')
instrument=$(awk -F, '/^instrument/ {print $2}' "$params_file" | tr -d '\r')
q_floor=$(awk -F, '/^q_floor/ {print $2}' "$params_file" | tr -d '\r')
q_cutoff=$(awk -F, '/^q_cutoff/ {print $2}' "$params_file" | tr -d '\r')
cutoff_pct=$(awk -F, '/^cutoff_pct/ {print $2}' "$params_file" | tr -d '\r')
sample_names_line=$(grep "sample_names" "$params_file" | tr -d '[:space:]')
reorg=$(awk -F, '/^reorganize/ {print $2}' "$params_file" | tr -d '\r')

cd ${path}/pipeline

# Run the executable with arguments
./merge_reads3.sh "${path}" "${pear}" "${pear_overlap}" "${pear_stattest}" "${pear_pvalue}" "${merge}" "${cpus}" "${memory}" "${instrument}" "${q_floor}" "${q_cutoff}" "${cutoff_pct}" "${sample_names_line}" "${reorg}"

# Create quality filtering C++ script
cat << EOF > QF3.cpp
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
using namespace std;

int q_floor = ${q_floor};			//q_floor set using params.csv
int q_cutoff = ${q_cutoff};			//q_cutoff set using params.csv
float cutoff_pct = ${cutoff_pct};	//cutoff_pct set using params.csv

int main()
{
	//Initialize variables and open file streams
	string inMPERs="combined.fastq";
	ifstream inRawReads(inMPERs.c_str());
	ofstream outGoodReads("good_reads.csv"), outPoorReads("poor_reads.csv");
	size_t split;				//position of the split between sequence & quality string
	double q, g, pct;			//Q-score, # bases > q_cutoff, % bases > q_cutoff
	int j;						//sequence position when looping through quality string
	string line, seq, qual;		//combined.fastq line, DNA sequence, quality string
	bool discarded;				//whether read contains bases < q_floor
	
	//Print quality filtering variables to the console
	cout << "           Quality filtering settings:" << '\n' << '\t';
	cout << "           Proportion >=Q" << q_cutoff << " must be above " << cutoff_pct << '\n' << '\t';
	cout << "           All bases must have a Q-score >= Q" << q_floor << '\n';
	
	while(getline(inRawReads,line))
	{
		split = line.find(" ");				//Get split position
		qual = line.substr(split + 1);		//Get quality sequence
		seq = line.substr(0, split);		//Get DNA sequence

		g = 0;
		discarded = false;
		
		for(j=0; j < qual.length(); j++)	//For each character in the quality string
		{
			q = double(qual[j]) - 33;		//Convert character to Q-score (ASCII - 33)
			if(q < q_floor)					//If base quality is < q_floor,
			{
				discarded = true;			//mark sequence as discarded
				break;						//and break the for loop (ignore further bases)
			} else
			{
				if(q >= q_cutoff)			//If base quality is >= q_cutoff
				{
					g = g + 1;				//Add one to g (good bases count)
				}
			}
		}
		if(discarded == false)				//If sequence isn't already marked discarded,
		{
			pct = g / qual.length();		//Calculate percent bases >= q_cutoff
			if(pct > cutoff_pct)			//If percent meets standards,
			{
				outGoodReads << seq << '\n';//Write seq to good_reads.csv
			} else							//If not,
			{								//Write seq to poor_reads.csv
				//Include "E2" classification and actual percent >=q_cutoff
				outPoorReads << seq << "," << "E2" << "," << pct << '\n';
			}
		} else								//If sequence is already marked discarded,
		{									//Write seq to poor_reads.csv
			//Include "E1" classification and position of first subthreshold base
			outPoorReads << seq << "," << "E1" << "," << j << '\n';
		}
	}
	inRawReads.close();
	outGoodReads.close();
	outPoorReads.close();
	
	cout << "           High quality reads written to 'good_reads.csv'" << '\n';
	cout << "           Low quality reads written to 'poor_reads.csv'" << '\n' << '\t';
	cout << "           'E1' - 1+ bases < Q" << q_floor << " (first base index shown)" << '\n' << '\t';
	cout << "           'E2' - < " << cutoff_pct << " Q" << q_cutoff << " (actual percentage shown)" << '\n';

	return 0;
}
EOF

# Compile the quality filtering script (static compiling allows this to run through a cluster job)
${compiler} -static-libstdc++ -o QF3.out QF3.cpp
#Thank Silas and Tony 
echo "Thanks Silas and Tony"
