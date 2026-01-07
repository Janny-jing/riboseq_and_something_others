#! /bin/sh

for i in *.fna;do
	n=${i/.fna/}
	tRNAscan-SE --max ${n}.fna  -o ${n}_result1.txt -f ${n}_result2.txt -m ${n}_result3.txt
done
