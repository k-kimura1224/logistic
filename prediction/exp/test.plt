set key off
set key outside
set style boxplot nooutliers
set terminal aqua
plot "./biodeg_para_1.bp" using (1.0):2:(0):1 with boxplot, \
     "./biodeg_para_2.bp" using (4.0):2:(0):1 with boxplot, \
     "./biodeg_para_3.bp" using (7.0):2:(0):1 with boxplot, \
