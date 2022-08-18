src_exts = .c .cpp
obj_ext = .obj
exe_ext = .exe
tool ?= MinGW
$(foreach ext,$(src_exts),$(eval src += $(wildcard *$(ext))))
obj := $(src:%=%$(obj_ext))
dep := $(obj:$(obj_ext)=.dep)
ifeq ($(tool),MinGW)
	CC = gcc -Wno-pragmas
	CXX = g++ -Wno-pragmas
	CPP = cpp
else ifeq ($(tool),LLVM)
	CC = clang --target=x86_64-windows-gnu
	CXX = clang++ --target=x86_64-windows-gnu
	CPP = clang --target=x86_64-windows-gnu
endif
proj = geometria_opengl
NULLSTR =
CXXFLAGS += -std=c++2a -pedantic -pedantic-errors -Wall -Wextra -g -ggdb
CXXFLAGS += -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization
CXXFLAGS += -Wformat=2 -Wmissing-declarations
CXXFLAGS += -Wmissing-include-dirs -Wold-style-cast
CXXFLAGS += -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion
CXXFLAGS += -Wsign-promo -Wstrict-overflow=5 -Wswitch-default
CXXFLAGS += -Wundef -Werror -Wno-unused-parameter -Wno-unused-variable

INCFLAGS += -I"D:/GLFW/include" -I"D:/glad/include" -I"D:/glm-0.9.9.3/glm"
LDFLAGS += -L"D:/GLFW/lib"
LDLIBS += -lglfw3 -lgdi32 -lopengl32

EXEC = $(proj)$(exe_ext)
OBJ_DIR = obj
obj := $(obj:%=$(OBJ_DIR)/%)
# TEMPS = $(dep)
DEP_DIR = deps
dep := $(dep:%=$(DEP_DIR)/%)

END = ** everything seems up-to-date

build: greetings $(EXEC) end

$(DEP_DIR)/%.dep: %
	@$(CPP) $< $(INCFLAGS) -MM -MT "$(@:$(DEP_DIR)/%.dep=$(OBJ_DIR)/%$(obj_ext))" >$@

-include $(dep)

$(EXEC): $(obj)
	@echo  -- Linking into $@
	@$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)
	@echo  ++ No errors
	$(eval END = ** All operations done)

$(OBJ_DIR)/%.cpp$(obj_ext): %.cpp
	@echo  -- Compiling $< with $(firstword $(CXX))
	@$(CXX) $(INCFLAGS) $(CXXFLAGS) -o $@ -c $<

$(OBJ_DIR)/%.c$(obj_ext): %.c
	@echo  -- Compiling $< with $(firstword $(CC))
	@$(CC) $(INCFLAGS) $(CFLAGS) -o $@ -c $<

.PHONY: greetings
greetings: 
	@echo  ** In project $(proj)...
	@mkdir $(OBJ_DIR) 2>nul ||:
	@mkdir $(DEP_DIR) 2>nul ||:

# .PHONY: clean
# clean:
# 	@del $(TEMPS) 2>nul

.PHONY: end
end:
	@echo  $(END)
