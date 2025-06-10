/*                                          
 * variant_file.cpp
 *
 *  Created on: Jan 28, 2019
 *      Authors: Andreea Dreau, https://github.com/adreau
 *               Frank Chan, https://github.com/evolgenomics
 */
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <utility>
#include <vector>
#include <gzstream.h>

using namespace std;

string barcode_A="../data/BC_A_H4.txt";   // ??? 7|8, 
string barcode_B="../data/BC_B.txt";      // ??? 6, tested against (6+1, 6)
string barcode_C="../data/BC_C_H4.txt";   // ??? 7
string barcode_D="../data/BC_D.txt";      // ??? 6, tested against (0,6)

string barcode_ME="../data/BC_ME.txt";    // 10
string barcode_PLATE="../data/Plate_BC_8.txt";  // 8

map<string,string> bc_A;
map<string,string> bc_B;
map<string,string> bc_C;
map<string,string> bc_D;
map<string,string> bc_PLATE;
map<string,string> bc_ME;

typedef pair<unsigned int, unsigned int> pair_int;


unsigned int edit_distance(const std::string& s1, const std::string& s2)
{
	const size_t len1 = s1.size(), len2 = s2.size();
	vector<std::vector<unsigned int>> d(len1 + 1, vector<unsigned int>(len2 + 1));

	d[0][0] = 0;
	for(unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
	for(unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

	for(unsigned int i = 1; i <= len1; ++i)
		for(unsigned int j = 1; j <= len2; ++j)
                      // note that std::min({arg1, arg2, arg3}) works only in C++11,
                      // for C++98 use std::min(std::min(arg1, arg2), arg3)
                      d[i][j] = min(min( d[i - 1][j] + 1, d[i][j - 1] + 1), d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) );
	return d[len1][len2];
}

string min_edit_distance(const std::string& s1, map<string,string> bc, string code_letter){

  int min=10000;
  string code_min_dist;
  string code_str;
  int ed,occ=0;
  map<string,string>::iterator it;
  for(it=bc.begin();it!=bc.end();it++){
    ed=edit_distance(s1,it->first);
    //cout<<it->first<<":"<<ed<<endl;
    if(min>ed){
      min=ed;
      code_min_dist=it->second;
      code_str=it->first;
      occ=1;
    }else
    if(min==ed){
      occ++;
    }
  }

	if(occ>1)
		code_min_dist=code_letter+"00";

  //cout<<"Occurences: "<<occ<<endl;
  //cout<<"Min dist: "<<min<<endl;
  //cout<<"Code:"<<s1<<" Corrected to:" <<code_str<<endl;
  return code_min_dist;

}
void load_barcodes(map<string,string> &bc_list,string file, int len){

    ifstream barcode_file(file.c_str());
    string line;
    while(getline(barcode_file, line)){
      bc_list.insert( std::pair<string,string> (line.substr(len+1),line.substr(0,len)) );
    }


}

void load_ME(map<string,string> &bc_list,string file){

    ifstream barcode_file(file.c_str());
    string line;
    while(getline(barcode_file, line)){
      bc_list.insert( std::pair<string,string> (line.substr(3),line.substr(0,2)) );
    }


}

void load_all_barcodes(){

  load_barcodes(bc_A,barcode_A, 3);
  load_barcodes(bc_B,barcode_B, 3);
  load_barcodes(bc_C,barcode_C, 3);
  load_barcodes(bc_D,barcode_D, 3);
  load_barcodes(bc_PLATE,barcode_PLATE, 4);
  load_ME(bc_ME,barcode_ME);

}

//getStagger(R2,staggerME,10, bc_ME);

void getStagger(string &R2_prefix, string &staggerME, 
						int code_total_length, map<string,string> bcM){

	string codeME_inFile;
  map<string,string>::iterator m;
 // for(int i=0;i<2;i++)
 //  getline(R2, line);

//	RX1=line;
  // ??? code_total_length == 10, matches BC_ME.txt
  if(R2_prefix.length()<code_total_length){
    staggerME="S9";   // ??? fails if #stagger entries > 8
  }else{
		codeME_inFile=R2_prefix.substr(0,10);
    m=bcM.find(codeME_inFile);
    if(m==bcM.end()){
      staggerME=min_edit_distance(codeME_inFile,bcM,"S");
    }else
      staggerME=m->second;
  }

	//getline(R2, line);
	//getline(R2, line);
	//QX1=line;
}

//getCode(I2,codeB,codeD,codeA,codeC,RX1,QX1,read_type1,read_type2,22,stagger,"B","D","A","C", bc_B, bc_D, bc_A, bc_C);
//getCode(I2,codeB,codeD,codeA,codeC,RX1,QX1,read_type1,read_type2,22,stagger,"B","D","A","C", bc_B, bc_D, bc_A, bc_C);

void getCode(igzstream &I2, string &R3_orig, string &R3_qual, string &codeB, string &codeD, string &codeA, string &codeC,
						string& RX1, string& QX1, string& read_type1, string& read_type2,
						int code_total_length, string &staggerME, string code_letter1, string code_letter2, string code_letter3, string code_letter4, 
						map<string,string> bc1, map<string,string> bc2, map<string,string> bc3, map<string,string> bc4){

  // ??? bc1 = bc_B, bc2 = bc_D, bc3 = bc_A, bc4 = bc_C
	read_type1="correct";
	read_type2="correct";   // ??? const 
  
  string line;
	string codeB_inFile,codeD_inFile,codeA_inFile,codeC_inFile;
  map<string,string>::iterator b;
  map<string,string>::iterator d;
  map<string,string>::iterator a;
  map<string,string>::iterator c;

  string stagger_passed;
  stagger_passed=staggerME.substr(1,1);
  stringstream sss; 
  int sstagger;	
  sss << stagger_passed;
	sss >> sstagger;
  for(int i=0;i<2;i++)
    getline(I2, line);
	line=line.append(R3_orig);
  if (sstagger > 3) {
    int dummy = 0;
  }

	RX1=line;

  if(line.length()<code_total_length){
    codeB=code_letter1+"00";
    codeD=code_letter2+"00";
    codeA=code_letter3+"00";
    codeC=code_letter4+"00";
  }else{
	//Deal with B first
    // ??? same in demult_fastq:
    // ??? one character skipped between codeB_inFile and codeD_inFile
    // ??? why? 
		codeB_inFile=line.substr(6+1,6);    // ??? unbounded in demult_fastq
    b=bc1.find(codeB_inFile);
    if(b==bc1.end()){
      codeB=min_edit_distance(codeB_inFile,bc1,code_letter1);
			read_type1="corrected";// ??? forced to "corrected", ed ignored, codeB used
    }else
      codeB=b->second;

	//Deal with D next
		codeD_inFile=line.substr(0,6);
    d=bc2.find(codeD_inFile);
    if(d==bc2.end()){
      codeD=min_edit_distance(codeD_inFile,bc2,code_letter2);
			read_type1="corrected";   // ??? forced to "corrected", ed ignored, codeD used
    }else
      codeD=d->second;

  // ??? codeB_inFile codeA_inFile flushed (nothing skipped).
	//Deal with the staggering_A
		codeA_inFile=line.substr(13,6+sstagger);  // ??? 13 = |bc1| + 1 + |bc2|
 //   cout<< "Passed stagger: " << stagger_passed << endl;
 //   cout<< "A match: " << line.substr(13,6+sstagger) << endl;
    assert(6+sstagger == codeA_inFile.length());
		a=bc3.find(codeA_inFile);
    if(a==bc3.end()){
      // ??? bc3 contains codes x of |6|, |7| and |8|
      // ??? tested against |6|
      // ??? ed >= |x| - |6|
      codeA=min_edit_distance(codeA_inFile,bc3,code_letter3);
			read_type1="corrected";   // ??? forced to "corrected", ed ignored, codeA used
    }else
      codeA=a->second;
	
  // ??? codeA_inFile codeC_inFile flushed (nothing skipped).
	//Deal with the staggered_C
		codeC_inFile=line.substr(19+sstagger+1,7);    // ??? missing in code_total_length
    //cout<< "C match: " << line.substr(19+sstagger+1) << endl;
    c=bc4.find(codeC_inFile);
    if(c==bc4.end()){
      codeC=min_edit_distance(codeC_inFile,bc4,code_letter4);
			read_type1="corrected";    // ??? forced to "corrected", ed ignored, codeC used
    }else
      codeC=c->second;

  }

	if(codeB.compare(code_letter1+"00")==0 || codeD.compare(code_letter2+"00")==0 || codeA.compare(code_letter3+"00")==0 || codeC.compare(code_letter4+"00")==0){
		read_type1="unclear";
	}

  getline(I2, line);
	getline(I2, line);
	QX1=line;
	QX1=QX1.append(R3_qual);
}

void getPlateCode(igzstream &I1, string &codePLATE,
                                                string& RX1, string& QX1, string& read_type3,
                                                int code_total_length, string code_letter5,
                                                map<string,string> bc5){

        read_type3="correct";

  string line;
        string codePLATE_inFile;
  map<string,string>::iterator plate;

  for(int i=0;i<2;i++)
    getline(I1, line);

        RX1=line;

  if(line.length()<code_total_length){
    codePLATE=code_letter5+"00";
  }else{

//NEW FC - I haven't changed this yet
        //Deal with PLATE aka E
      //FC - 20220324 - HARDCODE for 7 here as a HACK          
	  // - ORIGINAL - codePLATE_inFile=line.substr(0,8);
	  codePLATE_inFile=line.substr(0,code_total_length);
    plate=bc5.find(codePLATE_inFile);
    if(plate==bc5.end()){
      codePLATE=min_edit_distance(codePLATE_inFile,bc5,code_letter5);
                        read_type3="corrected";
    }else
      codePLATE=plate->second;

  }

        if(codePLATE.compare(code_letter5+"00") == 0){
                read_type3="unclear";
        }

  getline(I1, line);
        getline(I1, line);
        QX1=line;
}


int main (int argc, char* argv[])
{

  load_all_barcodes();
  cout << "loaded barcodes: " << bc_A.size() << " A, " << bc_B.size() << " B, "
                              << bc_C.size() << " C, " << bc_D.size() << " D, " << bc_PLATE.size() << " PLATE, " << bc_ME.size() << " ME " << endl;


  string path_to_reads= (argc < 2) ? "../data/" : argv[1];
  string path_output = (argc < 2) ? "../data/" : argv[2];
                            
  string R1_file=path_to_reads+"_gen_test_H4_R1_001.fastq.gz";
  string R2_file=path_to_reads+"_gen_test_H4_R4_001.fastq.gz";
  // string R3_file=path_to_reads+"test_H4_R3_001.fastq.gz";   // ??? doesn't exist
  string R3_file=path_to_reads+"_gen_test_R3_001.fastq.gz"; 
  string I1_file=path_to_reads+"_gen_test_H4_I1_001.fastq.gz";
  string I2_file=path_to_reads+"_gen_test_H4_R2_001.fastq.gz"; 

  string R1_outfile=path_output+"_gen__R1_001.fastq.gz";
  string R2_outfile=path_output+"_gen__R2_001.fastq.gz";

	string clearBC_logfile=path_output+"_clearBC.log";
	string unclearBC_logfile=path_output+"_unclearBC.log";

  igzstream R1(R1_file.c_str());
  igzstream R2(R2_file.c_str());
  igzstream R3(R3_file.c_str());
  igzstream I1(I1_file.c_str());
  igzstream I2(I2_file.c_str());

	ogzstream R1_out(R1_outfile.c_str());
	ogzstream R2_out(R2_outfile.c_str());


	ofstream clearBC_log(clearBC_logfile.c_str());
	ofstream unclearBC_log(unclearBC_logfile.c_str());


  string codeA,codeB,codeC,codeD,codePLATE;
  string R3_orig, R3_qual, R2_orig, R2_prefix, staggerME;
	string RX1, RX2, QX1, QX2, PX1, PX1Q;   // ??? RX2, QX2 not used
	//read type: correct, corrected, unclear
	string read_type1, read_type2, read_type3;

	string R1_out_name,R2_out_name;

	map<string,pair_int> clear_read_map;
	map<string,pair_int>::iterator it_clear;
	map<string,int> unclear_read_map;
	map<string,int>::iterator it_unclear;

  string line;
  string name;
  string code;
  
  int posName;
  stringstream ss; 
  int stagger=9;
  string stagger_num;
	while (getline(R1, line))
	{

	

    posName=line.find(" ");
  	name=line.substr(0,posName+1);
	
	//Just getting the header line from Read3 out of the way...
    getline(R3, R3_orig);

  // ??? Read2?
	//Get the first read from Read2 to figure out what the stagger should be
    getline(R3, R3_orig);
    getline(R3, R3_qual);
    getline(R3, R3_qual);

	//Just getting the header line from Read2 out of the way...
	getline(R2, line);
	//Get the first read from Read2 to figure out what the stagger should be
	getline(R2, R2_orig);
  // ??? 10 <-> BC_ME.txt code length
	R2_prefix=R2_orig.substr(0,10);   // ??? unchecked |R2_orig|
	getStagger(R2_prefix,staggerME,10, bc_ME);
    stagger_num=staggerME.substr(1,1);
	ss << stagger_num;
	ss >> stagger;


  // ??? hard coded #stagger
	if (stagger < 3) {
		getCode(I2, R3_orig, R3_qual, codeB,codeD,codeA,codeC,RX1,QX1,read_type1,read_type2,22,staggerME,"B","D","A","C", bc_B, bc_D, bc_A, bc_C);
//	    getCode(I2,codeB,codeD,RX2,QX2,read_type2,13, "B", "D", bc_B, bc_D);
	} else {
    // ??? never reached
		codeA="A00";
		codeB="B00";
		codeC="C00";
		codeD="D00";
		read_type1="unclear";
		read_type2="unclear";
		RX1="NNNNNNNNNNNNNNNNNNNNNNNNNNNN";   // ??? 28
		RX2="NNNNNNNNNNNNNNNNNNNNNNNNNNNN";   // ??? 28
		QX1="----------------------------";   // ??? 28
		QX2="----------------------------";   // ??? 28 unused
	}

		//append BX tag

       //20230324 - FC - hard code to change the length of BC_PLATE to 7
		//ORIGINAL - getPlateCode(I1,codePLATE,PX1,PX1Q,read_type1,8,"P",bc_PLATE);
       getPlateCode(I1,codePLATE,PX1,PX1Q,read_type3,8,"P",bc_PLATE);

    name=name.append("BX:Z:");
    name=name.append(codeA);
    name=name.append(codeC);
    name=name.append(codeB);
    name=name.append(codeD);
    name=name.append("-");
    name=name.append(codePLATE);

		code=codeA+codeC+codeB+codeD+codePLATE;

		//append RX tag
		name=name.append("\tRX:Z:");
		name=name.append(RX1);
		name=name.append("+");
                name=name.append(PX1);
		//name=name.append("+");
		//name=name.append(RX2);


		//append QX tag
		name=name.append("\tQX:Z:");
		name=name.append(QX1);
                name=name.append("+");
                name=name.append(PX1Q);
		//name=name.append("+");
		//name=name.append(QX2);
	
//		cout << endl;
//		cout << "Read2_prefix: " << endl;
//	cout << R2_prefix << endl;
//	cout << "Stagger: " << staggerME << " - " << stagger_num << endl;
//	cout << "Name: " << name << endl;

    R1_out<<name<<endl;
    R2_out<<name<<endl;
    //get the codeA value and determine clip size based on codeA integer
    string codeA_value = codeA.substr(1,3);
    unsigned int codeA_num = atoi(codeA_value.c_str());
    unsigned int clip_size=0;
    if (codeA_num > 0 &&  codeA_num <= 32)
       clip_size = 17;  // |A| 6
    else if (codeA_num > 32 && codeA_num <=64)
       clip_size = 18;  // |A| 7
    else
       clip_size = 19;  // |A| 0|8   

    string R2_orig_clipped = R2_orig.substr(clip_size, R2_orig.size()-clip_size+1);

        //Writing out the original Read2 without any modification;
    R2_out<<R2_orig_clipped<<endl;

    for(int i=0;i<3;i++){
      getline(R1, line);
      R1_out<<line<<endl;
    }
    for(int i=0;i<2;i++){
      getline(R2, line);
      if (i == 1) // clip the quality line identical to how it was done for R2_orig->R2_orig_clipped
          line = line.substr(clip_size, line.size()+1);
      R2_out<<line<<endl;
    }
		//sort reads clear vs unclear
    // ??? boring output, always empty (stagger < 3) or full of "A00B00C00D00..."
		if(read_type1.compare("unclear")==0 || read_type2.compare("unclear")==0){

			it_unclear=unclear_read_map.find(code);
			if(it_unclear!=unclear_read_map.end())
				it_unclear->second++;
			else
				unclear_read_map.insert( pair<string,int>(code,1));

    }else{

			it_clear=clear_read_map.find(code);
			if(read_type1.compare("corrected")==0 || read_type2.compare("corrected")==0){

				if(it_clear!=clear_read_map.end())
					it_clear->second.second++;
				else
					clear_read_map.insert( make_pair(code, make_pair(0,1)) );

	    } else{

				if(it_clear!=clear_read_map.end())
					it_clear->second.first++;
				else
					clear_read_map.insert( make_pair(code, make_pair(1,0)) );

			}
		}

  }

	R1_out.close();
  R2_out.close();


	clearBC_log << "Barcode \t Correct reads \t Corrected reads" <<endl;
	for ( map<string,pair_int>::iterator it=clear_read_map.begin(); it!=clear_read_map.end(); ++it)
		clearBC_log << it->first << "\t" << it->second.first << "\t" << it->second.second << endl;
	clearBC_log.close();


	unclearBC_log << "Barcode \t Reads" <<endl;
	for ( map<string,int>::iterator it=unclear_read_map.begin(); it!=unclear_read_map.end(); ++it)
		unclearBC_log << it->first << "\t" << it->second << endl;
	unclearBC_log.close();

  return 0;
}
