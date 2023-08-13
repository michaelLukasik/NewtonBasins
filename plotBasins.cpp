#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
int plotTest() {
    plt::plot({ 1,3,2,4 });
    plt::show();
}