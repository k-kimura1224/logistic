#!/bin/sh
rm -f *.log

rm -rf torima
mkdir torima

mv *.sol torima
mv torima/*init* .
rm -rf torima
