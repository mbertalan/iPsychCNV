#!/usr/bin/perl

print("Chr\tStart\tStop\tNumSNPs\tLength\tState\tCN\tFile\tstartsnp\tendsnp\tconf\n");
while(<STDIN>) # chr1:1219382-1219456          numsnp=3      length=75          state1,cn=0 /home/tspa/projects/cnvCall/inputCNV/9372161190_R04C02 startsnp=exm2062 endsnp=exm2067 conf=9.184
{	
	# chr15:42153994-42164078       numsnp=22     length=10,085      state5,cn=3 /home/tspa/projects/cnvCall/inputCNV/9448203044_R05C01 startsnp=exm1152505 endsnp=exm1152639 conf=21.722
	
	if($_ =~ /chr(\S+)\:(\d+)\-(\d+)\s+numsnp\=(\d+)\s+length\=(\S+)\s+(\S+),cn\=(\d+)\s+(\S+)\s+startsnp\=(\S+)\s+endsnp\=(\S+)\s+conf\=(\S+)/)
	{
		$Chr = $1;
		$Start = $2;
		$Stop = $3;
		$NumSNPs = $4;
		$length = $5;
		$State = $6;
		$CN = $7;
		$File = $8;
		$startsnp = $9;
		$endsnp = $10;
		$conf = $11;
		$length =~ s/,//g;
		print("$Chr\t$Start\t$Stop\t$NumSNPs\t$length\t$State\t$CN\t$File\t$startsnp\t$endsnp\t$conf\n");
	}
}
