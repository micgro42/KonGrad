#!/bin/bash
./makedoc.sh
chmod a+r ../doc/html/*

if [ "`hostname`" = "Leawyn" ]; then
scp -rpC ../doc/html/* micgro42@eiger.physik.hu-berlin.de:CP32014/KonGrad/trunk/doc/html
fi
# cp -rp ../doc/html/* $HOME/public_html/cp3/KonGrad/
