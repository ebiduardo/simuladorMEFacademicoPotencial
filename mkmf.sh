gnuMake() {
(cd bin; mkmf -t ../gnu-linux.mk -c"-DwithHYPRE" -p simulador.exe -v ../src)
mv bin/Makefile bin/MakefileG
(cd bin; make -f MakefileG clean)
(cd bin; make -f MakefileG)
mv bin/simulador.exe bin/simuladorG.exe
}

intelMake() {
(cd bin; mkmf -t ../intel-linux.mk -c"-DwithHYPRE" -p simulador.exe -v ../src)
mv bin/Makefile bin/MakefileI

(cd bin; make -f MakefileI clean)
(cd bin; make -f MakefileI)
mv bin/simulador.exe bin/simuladorI.exe
}

intelMake
#gnuMake
