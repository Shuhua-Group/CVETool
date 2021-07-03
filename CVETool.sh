#!/bin/bash

## Environments
bin="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source ${bin}/parameter.txt
WD=`pwd`

FLAG_bam_list=0
FLAG_root_list=0
FLAG_ref_file=0
FLAG_sv_region=0
FLAG_snp_pos=0
FLAG_exp_file=0
FLAG_cov_file=0
taskname="myfile"
taskcount=1

if [ x$1 != x ]
then
	while getopts "b:f:r:v:p:e:c:n:t:" arg
	do
		case $arg in
		b)
			bam_list=$OPTARG
			FLAG_bam_list=1
			;;
		r)
			root_list=$OPTARG
			FLAG_root_list=1
			;;
		f)
			ref_file=$OPTARG
			FLAG_ref_file=1
			;;
		v)
			sv_region=$OPTARG
			FLAG_sv_region=1
			;;
		p)
			snp_pos=$OPTARG
			FLAG_snp_pos=1
			;;
		e)
			exp_file=$OPTARG
			FLAG_exp_file=1
			;;
		c)
			cov_file=$OPTARG
			FLAG_cov_file=1
			;;
		n)
			taskname=$OPTARG
			;;
		t)
			taskcount=$OPTARG
			;;
		?)
			echo "FAIL: Unknown argument. Please check again."
			exit 1
			;;
		esac
	done

	#check parameters
	until [ $FLAG_bam_list -eq 1 -a $FLAG_ref_file -eq 1 -a $FLAG_sv_region -eq 1 -a $FLAG_snp_pos -eq 1 -a $FLAG_exp_file -eq 1 -a $FLAG_cov_file -eq 1 ]
	do
		echo "FAIL: Variables uncomplete. Please check again."
		exit 1
	done
	
	RA1="realpath "`echo $bam_list`
	bam_list=`$RA1`

	if [ $FLAG_root_list -eq 1 ]
	then
		RA2="realpath "`echo $root_list`
		root_list=`$RA2`
	fi

	RA3="realpath "`echo $exp_file`
	exp_file=`$RA3`

	RA4="realpath "`echo $cov_file`
	cov_file=`$RA4`

	#clean failed runs
	rm -r ${WD}/${taskname} 2>/dev/null
	mkdir ${WD}/${taskname}
	cd ${WD}/${taskname}
	WD=`pwd`
	
	sampleID=`cat ${bam_list} | cut -f 1`
	cat ${bam_list} | cut -f 1 > ${WD}/id.list
	tasks=0
	for id in ${sampleID[@]}
	do
		bamfile=`cat ${bam_list} | grep "^${id}	" | cut -f 2`
		if [ $FLAG_root_list -eq 1 ]
		then
			rootfile=`cat ${root_list} | grep "^${id}	" | cut -f 2`
			${python3} ${bin}/freebayes_cnvnator.py -d ${WD} -id ${id} -bam ${bamfile} -cr ${rootfile} -ref ${ref_file} -sv_region ${sv_region} -snp_region ${snp_pos} -freebayes_loc ${freebayes} -cnvnator_loc ${cnvnator} &>>${WD}/${id}.log &
		else
			${python3} ${bin}/freebayes_cnvnator.py -d ${WD} -id ${id} -bam ${bamfile} -ref ${ref_file} -sv_region ${sv_region} -snp_region ${snp_pos} -freebayes_loc ${freebayes} -cnvnator_loc ${cnvnator} &>>${WD}/${id}.log &
		fi
		echo "Processing ${id}..."
		tasks=`expr $tasks + 1`
		if [ $tasks -ge $taskcount ]
		then
			wait
			tasks=`expr $tasks - $taskcount`
		fi
	done
	wait

	${python3} ${bin}/freebayes_cnvnator_summary.py -d ${WD} -l ${WD}/id.list -o ${WD}

	${perl} ${bin}/regression.pl ${WD} snp_summary.txt sv_summary.txt ${exp_file} ${cov_file} ${bin}/regression.R ${Rscript}
	${perl} ${bin}/view.pl ${WD} SV_SNP.paired.tmp
	cat view.tmp0 | sed 's/	$//g' | sed 's/^	//g' > view.tmp
	${python3} ${bin}/make_view_file.py view.tmp view.txt
	${Rscript} ${bin}/view.R view.txt ${WD}/Rplots.pdf &>/dev/null

	cat sv_summary.txt | sed "s/chr	pos	end/#CHROM	POS	END/g" > ./../${taskname}.genotype.txt
	cat snp_summary.txt | sed -n "2p" | cut -f 1,2 > snp1.tmp
	echo "." > snp2.tmp
	cat snp_summary.txt | sed -n "2p" | cut -f 3- > snp3.tmp
	paste snp1.tmp snp2.tmp snp3.tmp >> ./../${taskname}.genotype.txt
	mv SV_SNP.paired.txt ./../${taskname}.stat.txt
	mv Rplots.pdf ./../${taskname}.Rplot.pdf
	rm *tmp *tmp0
	
	TIMENOW=`date`
	echo ${TIMENOW}"	Softname.sh complete!"
else
	echo "	CVETool.sh
	[USAGE]
	-b bam list file
	-r root list file(option)
	-r reference file
	-v SV region
	-p SNP position
	-e gene expression file
	-c covariates file
	-n taskname
	-t parallel task count
	example: CVETool.sh -b bam.list -r human_g1k_v37.fasta -v 3:81595569-81599665 -p 3:81597974 -e exp.txt -c cov.txt -n taskname -t 10"
fi
