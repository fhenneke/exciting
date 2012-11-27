


.NOTPARALLEL:

default: build/make.inc all 

all:    serial mpi smp spacegroup stateinfo stateconvert species

build/make.inc:
	perl ./setup.pl

include build/make.inc

serial:
	cd build/serial; $(MAKE) 

mpi:
	cd build/mpi; $(MAKE) 

smp:
	cd build/smp; $(MAKE)

debug:
	cd build/debug; $(MAKE)

mpiandsmp:
	cd build/mpiandsmp; $(MAKE)

test::
	cd test/; $(MAKE) summary



doc:  spacegroupdoc stateconvertdoc stateinfodoc inputdoc excitingfuncdoc Splitt_inputdoc speciesdoc 

excitingfuncdoc::
	$(MAKE) -f build/Make.common doc

spacegroupdoc::
	cd src/spacegroup; $(MAKE) doc;\
	mv spacegroup.pdf ../../docs/spacegroup
	xsltproc xml/schema/schemaexpand.xsl xml/schema/symmetries.xsd>  xml/sgroupinput.xsd
	cd docs/spacegroup; \
	xsltproc --stringparam importancelevels "spacegroup" ../../xml/schematolatex.xsl ../../xml/sgroupinput.xsd > spacegroupinput.tex; \
	xsltproc --stringparam importancelevels "spacegroup" ../../xml/schematowikidot.xsl ../../xml/sgroupinput.xsd > ../../xml/schema/wiki/spacegroup ;\
	pdflatex spacegroupinput.tex; \
	pdflatex spacegroupinput.tex; \

speciesdoc::
	cd docs/species;\
	xsltproc --stringparam importancelevels "spacegroup" ../../xml/schematolatex.xsl ../../xml/species.xsd > species.tex;\
	pdflatex species.tex;pdflatex species.tex;



expandedschema::
	xsltproc xml/schema/schemaexpand.xsl xml/schema/input.xsd >xml/excitinginput.xsd ;\

inputdoc::expandedschema
	cd docs/exciting/;\
	xsltproc --stringparam importancelevels "essential expert" ../../xml/schematolatex.xsl ../../xml/excitinginput.xsd >excitinginput.tex;\
	xsltproc --stringparam importancelevels "essential expert" ../../xml/schematowikidot.xsl ../../xml/excitinginput.xsd >inputref.wikidot;\
	pdflatex excitinginput.tex;\
	pdflatex excitinginput.tex

Splitt_inputdoc::
	cd xml/schema && $(MAKE)

inputdocwiki:xml/schema/*.xsd 
	cd xml/schema; $(MAKE) 


stateconvertdoc::
	cd src/stateconvert; $(MAKE) doc;\
	mv stateconvert.pdf ../../docs/stateconvert

stateinfodoc::
	cd src/stateinfo; $(MAKE) doc;\
	mv stateinfo.pdf ../../docs/stateinfo

eos::
	cd src/eos; $(MAKE)

spacegroup::
	cd src/spacegroup; $(MAKE)

stateinfo::
	cd src/stateinfo; $(MAKE)

stateconvert::
	cd src/stateconvert; $(MAKE)  

species::libs
	cd src/species; $(MAKE)

libs:
	cd build/serial; $(MAKE) libs

clean:
	cd build/serial; $(MAKE) clean cleanlibs
	cd build/mpi; $(MAKE) clean cleanlibs
	cd build/smp; $(MAKE) clean cleanlibs
	cd build/debug; $(MAKE) clean cleanlibs
	cd build/mpiandsmp; $(MAKE) clean cleanlibs
	cd test/build ;$(MAKE) clean
	cd src/eos; $(MAKE) clean
	cd src/spacegroup; $(MAKE) clean
	cd src/species; $(MAKE) clean
	cd src/src_vdwdf; $(MAKE) clean
	cd src/stateinfo; $(MAKE) clean
	cd src/stateconvert; $(MAKE) clean
	rm -f *.o *.mod *~ fort.* ifc* *.gcno *.exe exdg.*
	rm -f bin/*
	rm -f interfaces/*
	rm -f docs/exciting/*
	rm -f docs/spacegroup/*
	rm -f src/leblaiklib/*.o src/leblaiklib/*.a

libxcclean:
	cd src/libXC && make clean 

tgz::doc #libxcclean
	tar  --exclude-from=".gitignore"  --transform 's,^,/exciting/,' -c -v  -f ./exciting.tar  *
	tar    -r -v  -f ./exciting.tar  --transform 's,^,/exciting/,' .git/HEAD  .git/refs .git/packed-refs \
	test/test02/reference/ test/test08/reference/
	gzip  -f --best ./exciting.tar 
	du -h ./exciting.tar.gz 

tidy:
	perl setup.pl tidy $(MAKE) 

vdwdf:
	cd src/src_vdwdf
	$(MAKE)
