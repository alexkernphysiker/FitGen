add_library(FitGen SHARED abstract.cpp fit.cpp filter.cpp paramset.cpp initialconditions.cpp)
target_link_libraries(FitGen gpwrap)