OBJS = CDG.o \
       CFA.o \
       Field.o \
       Basis.o \
       Grid.o \
       Cell.o \
       Polygon.o \
       Triangle.o \
       Edge.o \

FLAGS = -g -Wall
CC = g++

cdgTest: CDGTest.o ${OBJS}
	${CC} -o cdgTest CDGTest.o ${OBJS} ${FLAGS}
CDGTest.o: CDG.cpp CDG.o CFA.o Field.o Grid.o Polygon.o Basis.o Triangle.o Edge.o
	${CC} -c CDGTest.cpp ${FLAGS}

cfaTest: CFATest.o ${OBJS}
	${CC} -o cfaTest CFATest.o ${OBJS} ${FLAGS}
CFATest.o: CFATest.cpp CFA.o Field.o Basis.o Grid.o Cell.o Polygon.o Triangle.o Edge.o
	${CC} -c CFATest.cpp ${FLAGS}

cdgSSTest: CDGSingleCellTest.o ${OBJS}
	${CC} -o cdgSSTest CDGSingleCellTest.o ${OBJS} ${FLAGS}
CDGSingleCellTest.o: CDGSingleCellTest.cpp CDG.o CFA.o Field.o Basis.o Grid.o Cell.o Polygon.o Triangle.o Edge.o
	${CC} -c CDGSingleCellTest.cpp ${FLAGS}

cfaSSTest: CFASingleCellTest.o ${OBJS}
	${CC} -o cfaSSTest CFASingleCellTest.o ${OBJS} ${FLAGS}
CFASingleCellTest.o: CFASingleCellTest.cpp CFA.o Field.o Basis.o Grid.o Cell.o Polygon.o Triangle.o Edge.o
	${CC} -c CFASingleCellTest.cpp ${FLAGS}

fieldTest: FieldTest.o ${OBJS}
	${CC} -o fieldTest FieldTest.o ${OBJS} ${FLAGS}
FieldTest.o: FieldTest.cpp Field.o Basis.o Grid.o Cell.o Polygon.o Triangle.o Edge.o
	${CC} -c FieldTest.cpp ${FLAGS}

quadTest: QuadTest.o ${OBJS}
	${CC} -o quadTest QuadTest.o ${OBJS} ${FLAGS}
QuadTest.o: QuadTest.cpp Field.o Basis.o Grid.o Cell.o Polygon.o Triangle.o Edge.o
	${CC} -c QuadTest.cpp ${FLAGS}

basisTest: BasisTest.o ${OBJS}
	${CC} -o basisTest BasisTest.o ${OBJS} ${FLAGS}
BasisTest.o: BasisTest.cpp CDG.o CFA.o Field.o Basis.o Grid.o Cell.o Polygon.o Triangle.o Edge.o
	${CC} -c BasisTest.cpp ${FLAGS}

gridTest: GridTest.o ${OBJS}
	${CC} -o gridTest GridTest.o ${OBJS} ${FLAGS}
GridTest.o: GridTest.cpp Grid.o Cell.o Polygon.o Triangle.o Edge.o
	${CC} -c GridTest.cpp ${FLAGS}

triangleTest: TriangleTest.o ${OBJS}
	${CC} -o triangleTest TriangleTest.o ${OBJS} ${FLAGS}
TriangleTest.o: TriangleTest.cpp Triangle.o Edge.o
	${CC} -c TriangleTest.cpp ${FLAGS}

CDG.o: CDG.cpp CDG.h CFA.o Field.o Basis.o Grid.o Cell.o Polygon.o Triangle.o Edge.o
	${CC} -c CDG.cpp ${FLAGS}
CFA.o: CFA.cpp CFA.h Field.o Basis.o Grid.o Cell.o Polygon.o Basis.o Triangle.o Edge.o
	${CC} -c CFA.cpp ${FLAGS}
Field.o: Field.cpp Field.h Basis.o Grid.o Cell.o Polygon.o Basis.o Triangle.o Edge.o
	${CC} -c Field.cpp ${FLAGS}
Basis.o: Basis.cpp Basis.h Cell.o Polygon.o Triangle.o Edge.o
	${CC} -c Basis.cpp ${FLAGS}
Grid.o: Grid.cpp Grid.h Cell.o Polygon.o Triangle.o Edge.o
	${CC} -c Grid.cpp ${FLAGS}
Cell.o: Cell.cpp Cell.h Polygon.o Triangle.o Edge.o
	${CC} -c Cell.cpp ${FLAGS}
Polygon.o: Polygon.cpp Polygon.h Triangle.o Edge.o
	${CC} -c Polygon.cpp ${FLAGS}
Triangle.o: Triangle.cpp Triangle.h Edge.o
	${CC} -c Triangle.cpp ${FLAGS}
Edge.o: Edge.cpp Edge.h
	${CC} -c Edge.cpp ${FLAGS}

clean:
	rm *.o *Test
