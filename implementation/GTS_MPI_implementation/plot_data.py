import numpy as np

from matplotlib import pyplot as plt


cores = [1,4,9,16]

times = [22.6,6.57,3.00,2.00]

ideal = [1,4,9,16]

timesnp = np.array(times)

speeds = 22.6/timesnp


# one core: 22.6s
# four core: 6.57s
# 9 core: 3.00s
# 16 core: 2.00s

#
# Plotting
#
# pd.read_csv(fpath,sep="\t") // DataFrame.from_csv

# f, (ax1, ax2, ax3) = plt.subplots(3, sharex=False, sharey=False)

# plt.subplots_adjust(hspace=.4)

# ax1.set_title("Initialization", loc="left")
# ax1.set_ylabel('Time (s)')

# ax1.errorbar(matrix_sizes, t_initialization_means, yerr=t_initialization_errors)

# ax2.set_title("Matrix Filling", loc="left")
# ax2.set_ylabel('Time (s)')

# ax2.errorbar(matrix_sizes, t_matrixfill_means, yerr=t_matrixfill_errors)

# ax3.set_title("DGEMM", loc="left")
# ax3.set_ylabel('Time (s)')
# ax3.set_xlabel('Matrix Size', fontsize=12)


# ax3.errorbar(matrix_sizes, t_gemm_means, yerr=t_gemm_errors)
# plt.savefig(plot_path)

f, ax = plt.subplots(1, sharex=False, sharey=False)

plt.ylim((0,16))
plt.xlim((0,16))

ax.plot(cores,speeds ,linewidth=2, color="b")
ax.plot(cores,ideal,color="g",linestyle="dashed",linewidth=2)

ax.legend(loc="upper left", fontsize=10, labels=['Actual performance','Ideal performance'],prop={'size':15})

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(16) 
    # specify integer or one of preset strings, e.g.
    #tick.label.set_fontsize('x-small') 
    #tick.label.set_rotation('vertical')

for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(16) 

#ax.set_title("DGEMM", loc="left")
ax.set_ylabel('Performance (speedup)',fontsize=16)
ax.set_xlabel('Number of cores', fontsize=16)

plt.show()

plt.close()