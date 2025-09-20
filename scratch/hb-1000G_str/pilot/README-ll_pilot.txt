I have the folders as follows:
- SCRIPTS - central location for running SCRIPTS
- out - outputs
- job_out - output log from the jobs

14/05/2024 - FC
This is the folder where I'm processing the LL haplotag pilot dataset.
I'm first running a demultpliexing step to extract the molecular barcode. 

Getting everything to run on habrok isn't the most obvious - but it's a good learning experience anyway. Plus I'll need to do this at some point.

Other folders will be generated on the fly - and I'll try to document them.

#18.May.2024 - split files for easy demulitplexing:
for i in {1..4}; do zcat Pilot-1_R${i}_001.fastq.gz | split -l 467250466 --numeric-suffix - --filter='gzip > $FILE.gz' split_out/Pilot_1_R${i}_001.part_ --suffix-length=3 --additional-suffix=.fastq; done

21/May/2024 - FC
This is super stupid. I was trying to get things going on the job servers and it was not generating any results, so I was testing out the commands with a tmp/ folder.

Then I managed to unlink everything because I typed:
rm tmp/ * 
extra space. That's that. Grrrr.....

Quick - going back onto the Azenta website to download the original files. I'll need to re-do the 4 segments of the data for the demultiplexing.

Also, I need to grab the barcodes again.


An old file list contained:
[p314775@login1 pilot]$ ls -ltr
total 2402341012
-rwxr-----+ 1 p314775 hb-1000G_str           66 Apr 25 15:40 LL-GoNL-Pilot-1_R1_001.fastq.gz.md5
-rwxr-----+ 1 p314775 hb-1000G_str 528519830590 Apr 25 15:40 LL-GoNL-Pilot-1_R1_001.fastq.gz
-rwxr-----+ 1 p314775 hb-1000G_str           66 Apr 25 22:20 LL-GoNL-Pilot-1_R2_001.fastq.gz.md5
-rwxr-----+ 1 p314775 hb-1000G_str 606082691344 Apr 25 22:20 LL-GoNL-Pilot-1_R2_001.fastq.gz
-rwxr-----+ 1 p314775 hb-1000G_str           66 Apr 30 14:42 LL-GoNL-Pilot-1_I1_001.fastq.gz.md5
-rwxr-----+ 1 p314775 hb-1000G_str           66 Apr 30 14:42 LL-GoNL-Pilot-1_I2_001.fastq.gz.md5
-rwxr-----+ 1 p314775 hb-1000G_str  39100366402 May  1 14:04 LL-GoNL-Pilot-1_I2_001.fastq.gz
-rwxr-----+ 1 p314775 hb-1000G_str  83371694402 May 13 11:34 LL-GoNL-Pilot-1_I1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str        38308 May 14 10:39 test_R2_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str        24048 May 14 10:39 reads.grep
-rw-rw----+ 1 p314775 hb-1000G_str        20040 May 14 10:41 reads.grep.1
-rw-rw----+ 1 p314775 hb-1000G_str        34273 May 14 10:41 test_R1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str         6020 May 14 10:41 test_I1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str         2654 May 14 10:41 test_I2_001.fastq.gz
-rwxr-x---+ 1 p314775 hb-1000G_str           42 May 14 11:27 BC_ME.txt
-rwxr-x---+ 1 p314775 hb-1000G_str         1056 May 14 11:27 BC_D.txt
-rwxr-x---+ 1 p314775 hb-1000G_str         1152 May 14 11:27 BC_C_H4.txt
-rwxr-x---+ 1 p314775 hb-1000G_str         1056 May 14 11:27 BC_B.txt
-rwxr-x---+ 1 p314775 hb-1000G_str         1152 May 14 11:27 BC_A_H4.txt
-rwxr-x---+ 1 p314775 hb-1000G_str          336 May 14 11:29 Plate_BC_8.txt
-rwxr-x---+ 1 p314775 hb-1000G_str          156 May 14 11:29 Plate_BC_7.txt
-rw-rw----+ 1 p314775 hb-1000G_str        38308 May 14 11:32 test_R3_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str        36430 May 14 12:03 test_H4_R4_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str         8605 May 14 12:06 test_H4_R2_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str        34273 May 14 12:06 test_H4_R1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str         6020 May 14 12:06 test_H4_I1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str        36430 May 14 12:12 Haplo4_R3_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str        34273 May 14 12:15 Haplo4_R1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str         6020 May 14 12:15 Haplo4_I1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str         8605 May 14 12:16 Haplo4_R2_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str          216 May 14 12:24 Haplo4_focal_R1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str          144 May 14 12:27 Haplo4_focal_I1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str          179 May 14 12:29 Haplo4_focal_R2_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str          157 May 14 12:30 Haplo4_focal_R4_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str          151 May 14 12:30 Haplo4_focal_R3_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str           16 May 14 12:30 out_unclearBC.log
-rw-rw----+ 1 p314775 hb-1000G_str          197 May 14 12:30 out_R2_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str          205 May 14 12:30 out_R1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str           63 May 14 12:30 out_clearBC.log
-rw-rw----+ 1 p314775 hb-1000G_str        34273 May 14 14:50 Haplo4_small_R1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str         2679 May 14 14:50 Haplo4_small_I1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str        36430 May 14 14:50 Haplo4_small_R4_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str         6159 May 14 14:50 Haplo4_small_R3_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str        11665 May 14 14:50 Haplo4_small_R2_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str          828 May 14 14:52 small_unclearBC.log
-rw-rw----+ 1 p314775 hb-1000G_str        45424 May 14 14:52 small_R2_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str        49583 May 14 14:52 small_R1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str         9512 May 14 14:52 small_clearBC.log
-rw-rw----+ 1 p314775 hb-1000G_str  35884851935 May 14 18:05 Pilot-1_I1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str            0 May 14 22:30 Pilot-1_I1_001.fastq.linecount
-rw-rw----+ 1 p314775 hb-1000G_str            0 May 14 22:30 LL-GoNL-Pilot-1_I2_001.fastq.linecount
lrwxrwxrwx. 1 p314775 hb-1000G_str           59 May 16 10:39 Pilot-1_R1_001.fastq.gz -> /scratch/hb-1000G_str/pilot/LL-GoNL-Pilot-1_R1_001.fastq.gz
drwxrws---+ 2 p314775 hb-1000G_str         4096 May 16 11:57 out
-rwxrwxrwx. 1 p314765 p314765       80017348852 May 16 15:34 Pilot-1_R3_001.fastq.gz
-rwxrwxrwx. 1 p314765 p314765      152375489809 May 17 00:01 Pilot-1_R2_001.fastq.gz
-rw-rw-r--. 1 p314765 p314765      524807037662 May 17 22:50 Pilot-1_R4_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str          779 May 18 22:23 README-ll_pilot.txt
-rw-rw----+ 1 p314775 hb-1000G_str            0 May 19 23:36 LL-GoNL-Pilot-1_demux_unclearBC.log
-rw-rw----+ 1 p314775 hb-1000G_str            0 May 19 23:36 LL-GoNL-Pilot-1_demux_clearBC.log
-rw-rw----+ 1 p314775 hb-1000G_str            0 May 19 23:44 LL-GoNL-Pilot-1_demult_unclearBC.log
-rw-rw----+ 1 p314775 hb-1000G_str            0 May 19 23:44 LL-GoNL-Pilot-1_demult_clearBC.log
drwxrws---+ 2 p314775 hb-1000G_str         4096 May 21 11:59 job_out
drwxrws---+ 2 p314775 hb-1000G_str         4096 May 21 13:21 split_out
drwxrws---+ 3 p314775 hb-1000G_str         4096 May 21 14:00 SCRIPTS
drwxrws---+ 2 p314775 hb-1000G_str         4096 May 21 14:00 tmp
-rw-rw----+ 1 p314775 hb-1000G_str  98012160000 May 21 14:03 LL-GoNL-Pilot-1_demux_R2_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str 107020271616 May 21 14:03 LL-GoNL-Pilot-1_demux_R1_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str  97911660544 May 21 14:03 LL-GoNL-Pilot-1_demult_R2_001.fastq.gz
-rw-rw----+ 1 p314775 hb-1000G_str 106909024256 May 21 14:03 LL-GoNL-Pilot-1_demult_R1_001.fastq.gz


21-05-2024/FC:
This is what I got from Marek regarding the reads and the barcode arrangement. *VERY* confusing!!!!! Note to self - need to have some staff to help me sort this out for real going forward.

--use-bases-mask=Y151,y13,i8y16,Y145
Haplo4_R1 is R1
HAplo4_R2 is the 13bp I1
Haplo4_R3 is the 16bp from the beginning of R2
Haplo4_R4 is the R2 minus 16bp
and the Haplo4_I1 is actually the I2

So:
ln -s R1 Haplo4_R1;
cutadapt -l 13 I1 Haplo4_R2
cutadapt -l 16 R2 Haplo4_R3
cutadapt -u 16 R2 Haplo4_R4


22-05-2024/FC: 
OK, so it wasn't all that bad. Overnight I was able to re-download the data from Azenta. Also through various transferring efforts [say, to MPI Tü], I was able to salvage some additional large files as well.

I have also established now that my demultiplexing code is working correctly, except that my earlier combination of I1+R2 wasn't quite right. Luckily, it was just redundant information, rather than wrong info.

So taking all of those together, this is some progress.

I'm now splitting the files into 60 chunks. It was a bit mysterious to me earlier why 60 chunks would be an even divider, until I realised that the flow cell has 120 tiles. So that makes sense now.

OK - so as of 8:42am, I have submitted the splitting jobs for R3 and R4. Once these are completed, I should be able to run the demultiplexing jobs for all the reads, split into 60 chunks, so that should be much faster.

Note also to self and GL - the I/O limit for the login node is IMMENSE. So yes even though it sounds silly to always copy the whole working file across to the node, it works out to be faster after all.Live and learn, I suppose.

24-05-2024/FC - 9am:
State-of-play: OK - so I did another dummy mistake. I was spitting files into 60 chunks - BUT - I forgot that the line count has to be divisible by 4... so. DUH! So now we end up having the situation where my demult script throws an error towards the end of the file and dumped the memory core. Arrrrggghhh.

So I re-started the split script yesterday at 3:20pm. It's still running - looks like another 4 hours or so? So it's almost 24h for the split script to run. Yikes.

In the meanwhile I built up the script for the next step, to split the files into each individual sample for easier processing down the line. I'm testing it on some of the first split-and-demult'd files just to see.

To run:
sbatch SCRIPTS/demult.sh Pilot_1 Pilot_1_BX

sbatch SCRIPTS/split_id.sh Pilot_1_BX Pilot_1_BX

#25-05-2024/FC - 11pm:
It seems to all run fine now. I'm going to have them grouped back together to make a single per-sample file.

#First step for gather:
cd /scratch/hb-1000G_str/pilot/split_demult/perSample
for i in {00..96}; do for n in N501 N502 P999 P000; do mkdir C${i}-$n; done; done
for i in {00..96}; do for n in N501 N502 P999 P000; do mv ../per_sample/*C${i}-$n* C${i}-$n/; done; done
mv perSample per_sample


#Generated per sample fastq file list like so:
for read in 1 2; do for c in {01..96}; do for p in N501 N502; do echo split_demult/per_sample/C${c}-$p/Pilot_1_BX_R${read}_001.C${c}-$p.fastq.gz; done; done; done > per_sample.fastq.filelist
#The full list should be like:
for read in 1 2; do for c in {00..96}; do for p in N501 N502 P000 P999; do echo split_demult/per_sample/C${c}-$p/Pilot_1_BX_R${read}_001.C${c}-$p.fastq.gz; done; done; done > per_sample.fastq.filelist.full

#looks like the collecting by individual went quite well.
#---- NO NO ---- there were some OOM kills it turns out. So I'll re-run the whole thing.
for i in `ls split_demult/per_sample | grep -v "P000"`; do sbatch SCRIPTS/collect.sh $i; done

#So now I'm trying to do cutadapt on the reads, so that we trim away all the adaptor sequences as well. This helps to avoid any mapping artefact due to adapter sequences.
sbatch SCRIPTS/cutadapt.arrayjob.Haplo4.sh



ls fastq | grep -v P000  | grep -v 000 | grep -v 999 > sample_ID.list

#30-May-2024/FC  - Checked to make sure all the input files are EMA-placed. It looks good!

[p314775@login1 bamfiles]$ ls  Pilot_1_BX.* | sed 's/\./\t/g' | grep BXnum | datamash groupby 3 count 1  
C01-N501	50
C01-N502	50
C02-N501	50
C02-N502	50
C03-N501	50
C03-N502	50
C04-N501	50
C04-N502	50
C05-N501	50
C05-N502	50
C06-N501	50
C06-N502	50
C07-N501	50
C07-N502	50
C08-N501	50
C08-N502	50
C09-N501	50
C09-N502	50
C11-N501	50
C11-N502	50
C12-N501	50
C12-N502	50
C13-N501	50
C13-N502	50
C14-N501	50
C14-N502	50
C15-N501	50
C15-N502	50
C16-N501	50
C16-N502	50
C17-N501	50
C17-N502	50
C18-N501	50
C18-N502	50
C19-N501	50
C19-N502	50
C20-N501	50
C20-N502	50
C21-N501	50
C21-N502	50
C22-N501	50
C22-N502	50
C23-N501	50
C23-N502	50
C24-N501	50
C24-N502	50
C25-N501	50
C25-N502	50
C26-N501	50
C26-N502	50
C27-N501	50
C28-N501	50
C28-N502	50
C29-N501	50
C29-N502	50
C30-N501	50
C30-N502	50
C31-N501	50
C32-N501	50
C32-N502	50
C33-N501	50
C33-N502	50
C34-N501	50
C34-N502	50
C35-N501	50
C35-N502	50
C36-N501	50
C36-N502	50
C37-N501	50
C37-N502	50
C38-N501	50
C38-N502	50
C39-N501	50
C39-N502	50
C40-N501	50
C40-N502	50
C41-N501	50
C41-N502	50
C42-N501	50
C42-N502	50
C43-N501	50
C43-N502	50
C44-N501	50
C44-N502	50
C45-N501	50
C45-N502	50
C46-N501	50
C46-N502	50
C47-N501	50
C47-N502	50
C48-N501	50
C48-N502	50
C49-N501	50
C49-N502	50
C50-N501	50
C50-N502	50
C51-N501	50
C51-N502	50
C52-N501	50
C52-N502	50
C53-N501	50
C53-N502	50
C54-N501	50
C54-N502	50
C55-N501	50
C55-N502	50
C56-N501	50
C56-N502	50
C57-N501	50
C57-N502	50
C58-N501	50
C58-N502	50
C59-N501	50
C59-N502	50
C60-N501	50
C60-N502	50
C61-N501	50
C61-N502	50
C62-N501	50
C62-N502	50
C63-N501	50
C63-N502	50
C64-N501	50
C64-N502	50
C65-N501	50
C65-N502	50
C66-N501	50
C66-N502	50
C67-N501	50
C67-N502	50
C68-N501	50
C68-N502	50
C69-N501	50
C69-N502	50
C70-N501	50
C70-N502	50
C71-N501	50
C71-N502	50
C72-N501	50
C72-N502	50
C73-N501	50
C73-N502	50
C74-N501	50
C74-N502	50
C75-N501	50
C75-N502	50
C76-N501	50
C76-N502	50
C77-N501	50
C77-N502	50
C78-N501	50
C78-N502	50
C79-N501	50
C79-N502	50
C80-N502	50
C81-N501	50
C81-N502	50
C82-N501	50
C82-N502	50
C83-N501	50
C83-N502	50
C84-N501	50
C84-N502	50
C85-N501	50
C85-N502	50
C86-N501	50
C86-N502	50
C87-N501	50
C87-N502	50
C88-N501	50
C88-N502	50
C89-N501	50
C89-N502	50
C90-N501	50
C90-N502	50
C91-N501	50
C91-N502	50
C92-N501	50
C92-N502	50
C93-N501	50
C93-N502	50
C94-N501	50
C94-N502	50
C95-N501	50
C95-N502	50
C96-N501	50
C96-N502	50



