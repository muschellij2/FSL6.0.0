

# A Makefile for the miscmaths class

include ${FSLCONFDIR}/default.mk

PROJNAME = ciftiio

USRINCFLAGS = -g -DCIFTILIB_USE_XMLPP -I${FSLEXTINC} -I${INC_XML2} -I${INC_XML++} -I${INC_XML++CONF} -I${INC_BOOST} -I./Common -I./Nifti -I./Cifti 
USRLDFLAGS =  -L${LIB_NEWMAT}
COMMONOBJS = Common/AString.o Common/BinaryFile.o Common/MathFunctions.o Common/XmlAdapter.o Common/CiftiException.o Common/FloatMatrix.o Common/Vector3D.o	
NIFTIOBJS = Nifti/NiftiHeader.o Nifti/NiftiIO.o
CIFTIOBJS = Cifti/CiftiBrainModelsMap.o Cifti/CiftiFile.o Cifti/CiftiLabelsMap.o Cifti/CiftiMappingType.o Cifti/CiftiParcelsMap.o Cifti/CiftiScalarsMap.o Cifti/CiftiSeriesMap.o Cifti/CiftiVersion.o Cifti/CiftiXML.o Cifti/Label.o Cifti/LabelTable.o Cifti/MetaData.o Cifti/StructureEnum.o Cifti/VolumeSpace.o

OBJS = ${COMMONOBJS} ${NIFTIOBJS} ${CIFTIOBJS}
LIBS = -lcifti -lxml++-2.6 -lxml2 -lboost_filesystem -lboost_system -liconv -lz -lm 

%.o: %.cxx 
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

HFILES = */*.h
# The target "all" should ALWAYS be provided
#   typically it will just be another target name

all: libcifti.a

libcifti.a: ${OBJS}
	${AR} -r libcifti.a ${OBJS}

test : libcifti.a example/rewrite.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ example/rewrite.o ${LIBS}

clean:
	${RM} -f /tmp/fslgrot ${OBJS} example/*.o *.a *.exe core depend.mk