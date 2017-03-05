
CXX=g++

SRC_DIR=source/
HEADER_DIR=include/
OBJ_DIR=obj/
BIN_DIR=bin/

TARGET=$(BIN_DIR)MeshViewer

SRC=$(wildcard ./*.c ./*.cpp $(SRC_DIR)*.c $(SRC_DIR)*.cpp)
HEADER=$(wildcard ./*.h $(HEADER_DIR)*.h)
OBJ=$(addprefix $(OBJ_DIR),$(patsubst %.cpp,%.o,$(notdir $(SRC))))
#OBJ=$(patsubst %.cpp,%.o,$(notdir $(SRC)))

INCLUDE=$(addprefix -I,$(HEADER_DIR))

# opengl params
GL_LIBS= -I/usr/include -L/usr/local/lib -L/usr/X11R6/lib
GL_FLAGS= -lglut -lGLU -lGL -lX11 -lXext  -lXi -lm

$(TARGET):$(OBJ)
	@mkdir $(BIN_DIR) -p
	$(CXX) $(GL_LIBS) $^ -o $@ $(GL_FLAGS)
	@echo done
$(OBJ):$(OBJ_DIR)%.o:$(SRC_DIR)%.cpp
	@mkdir $(OBJ_DIR) -p
	$(CXX) -c $< -o $@ $(HEADER_DIR)
.PHONY:clean
clean:
	rm -rf $(TARGET) $(OBJ)

