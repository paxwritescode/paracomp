#ifndef STUDY_H
#define STUDY_H

#include <stdio.h>
#include <omp.h>

#include "generate.h"
#include "picard.h"
#include "matrix_tools.h"

void study_time_vs_size(void);
void study_dependencies_on_threads(void);

#endif