#!/bin/bash

for file in pdf/*.pdf
do
    base=$(basename $file .pdf)
    pdf2swf $file -o swf/$base.swf
    swfcombine -X 600 -Y 800 swf/rfxview.swf swf=swf/$base.swf -o swf/$base.swf
done


