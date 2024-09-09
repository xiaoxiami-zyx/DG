# Makefile template
# Compiler
CXX = icpx
# Compiler flags
CXXFLAGS = -Wall -Wextra -O3 -fPIC -xHost -ftree-vectorize -mavx -msse2 -std=c++23
# Dependency flags
DEPFLAGS = -MMD -MP

# openMP
# OPENMP = -qopenmp

# Source files
SRCS = mesh.cpp main.cpp param.cpp field.cpp gauss_Lobatto.cpp solver.cpp base_function.cpp \
		output.cpp
# Object files
OBJS = $(SRCS:.cpp=.o)
# Dependency files
DEPS = $(SRCS:.cpp=.d)

# 3rd
# gmsh_DIR = /home/xxm/gmsh
# eigen_Dir = /home/xxm/eigen
# amgcl_Dir = /home/xxm/amgcl

# Include paths
# LIB = -L$(gmsh_DIR)/lib
# INCLUDE = -I$(gmsh_DIR)/include -I$(eigen_Dir)/Eigen -I$(amgcl_Dir)

# LINK = -lgmsh

# Executable
TARGET = myprogram
# 生成可执行文件
$(TARGET): $(OBJS)
	$(CXX) $(OPENMP) $(CXXFLAGS) $(LIB) $^ -o $@ $(LINK)

# 生成目标文件
%.o: %.cpp
	$(CXX) $(OPENMP) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDE) -c $< -o $@

# 清理
clean:
	rm -f $(OBJS) $(TARGET) $(DEPS)

# 包含所有的依赖文件
-include $(DEPS)