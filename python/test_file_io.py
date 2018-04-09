import matplotlib.pyplot as plt

# import the package locally, support python 2 and 3 using if and else
if __name__ == '__main__':
    from read_pillar_bin import PillarTrackerBinReader
else:
    from .read_pillar_bin import PillarTrackerBinReader

# specify a path to a pillar binary file
fpath = r"C:\Users\xiaochun\Desktop\MDA_pillar_10sec_15min__R3D.bin"

# plot the first frame
reader = PillarTrackerBinReader()
[meta, data] = reader.loadfile(fpath)
plt.scatter(data.trackX[0, :], data.trackY[0, :])
plt.show()