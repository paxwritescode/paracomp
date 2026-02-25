 #ifndef STUDY_H
 #define STUDY_H

#include "seidel.h"
#include "def.h"

void study_dependencies_on_nodes(void);
void study_dependencies_on_threads(int rank, int size);

 #endif