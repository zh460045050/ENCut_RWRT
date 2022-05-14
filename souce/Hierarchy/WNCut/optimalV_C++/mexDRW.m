mex -c ../optimalB_C++/optimal_b.cpp optimal_v.cpp -I'/usr/local/include' -I'../optimalB_C++' -L'/usr/local/lib' -lopencv_core.2.4.13 -lopencv_highgui.2.4.13 -lopencv_imgproc.2.4.13

mex -O mex_optimalVH.cpp -I'/usr/local/include' -I'../optimalB_C++' -L'/usr/local/lib' -lopencv_core.2.4.13 -lopencv_highgui.2.4.13 -lopencv_imgproc.2.4.13  optimal_v.o optimal_b.o