#ifndef UTIL_HPP
#define UTIL_HPP

#if NDEBUG
#define DEBUG_PRINT(fmt, args...) /* Don't do anything in release builds */
#define DEBUG_PRINTLN(fmt, args...) /* Don't do anything in release builds */
#else
#define DEBUG_PRINT(fmt, args...) fprintf(stderr, "DEBUG: %s:%d:%s(): " fmt, \
__FILE__, __LINE__, __func__, ##args)
#define DEBUG_PRINTLN(fmt, args...) fprintf(stderr, "DEBUG: %s:%d:%s(): \n" fmt, \
__FILE__, __LINE__, __func__, ##args)
#endif

#include <battery/tuple.hpp>
#include <battery/vector.hpp>

template<class U>
void print_matrix(battery::vector<battery::vector<U> > matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        printf("[");
        for (int j = 0; j < matrix[i].size(); j++) {
            matrix[i][j].print();
            if (j != matrix[i].size() - 1) {
                printf(", ");
            }
        }
        printf("]\n");
    }
}


#endif //UTIL_HPP
