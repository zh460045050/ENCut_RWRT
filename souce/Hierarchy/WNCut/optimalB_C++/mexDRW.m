mex -c optimal_b.cpp -I'/usr/local/include' -L'/usr/local/lib' -lopencv_core.2.4.13 -lopencv_highgui.2.4.13 -lopencv_imgproc.2.4.13

mex -O mex_optimalb.cpp -I'/usr/local/include' -L'/usr/local/lib' -lopencv_core.2.4.13 -lopencv_highgui.2.4.13 -lopencv_imgproc.2.4.13  optimal_b.o