import numpy as np


class SteppingData(object):
    def __init__(self, dataFile):
        self.dataFile = dataFile
        self.rawData = np.loadtxt(self.dataFile)
        # assert(np.shape(self.rawData)[1] == 4)  # guarantee 4 columns in data file
        self.bindTimes = self.rawData[:, 1]
        self.unbindTimes = self.rawData[:, 0]
        self.nbx_bind = np.around(self.rawData[:, 2], decimals=12)
        self.fbx_bind = np.around(self.rawData[:, 3], decimals=12)
        self.nmx_unbind = np.around(self.rawData[:, 4], decimals=12)
        self.fmx_unbind = np.around(self.rawData[:, 5], decimals=12)
        self.nmx_bind = np.around(self.rawData[:, 6], decimals=12)
        self.fmx_bind = np.around(self.rawData[:, 7], decimals=12)

        self.near_step_len = self.nbx_bind[1:]-self.nbx_bind[:-1]
        self.far_step_len = self.fbx_bind[1:]-self.fbx_bind[:-1]

        self.onebound_times = self.bindTimes[1:]-self.unbindTimes[1:]
        self.bothbound_times = self.unbindTimes[1:]-self.bindTimes[:-1]

        self.initial_displacements = []
        self.final_displacements = []

        assert(len(self.nbx_bind)==len(self.fbx_bind))
        for s in range(1, len(self.nbx_bind)):
            # print('near', self.nbx_bind[s-1], self.nbx_bind[s],  self.nbx_bind[s-1]-self.nbx_bind[s])
            # print('far', '%.19g' % self.fbx_bind[s-1], '%.19g' % self.fbx_bind[s],
            #       self.fbx_bind[s-1]-self.fbx_bind[s])
            assert((self.nbx_bind[s-1] == self.nbx_bind[s])
                   or (self.fbx_bind[s-1] == self.fbx_bind[s]))
            if(self.nbx_bind[s-1]== self.nbx_bind[s] and self.fbx_bind[s-1] == self.nbx_bind[s]):
                continue
            if not (self.nbx_bind[s-1] == self.nbx_bind[s]):
                self.initial_displacements.append(self.nbx_bind[s-1]-self.fbx_bind[s-1])
                self.final_displacements.append(self.nbx_bind[s]-self.fbx_bind[s])
            elif not (self.fbx_bind[s-1] == self.fbx_bind[s]):
                self.initial_displacements.append(self.fbx_bind[s-1]-self.nbx_bind[s-1])
                self.final_displacements.append(self.fbx_bind[s] - self.nbx_bind[s])
        self.initial_displacements = np.asarray(self.initial_displacements)
        self.final_displacements = np.asarray(self.final_displacements)

        self.step_lengths = self.final_displacements - self.initial_displacements
        self.step_times = self.onebound_times + self.bothbound_times


class MovieData: 
    pass


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        data = SteppingData(sys.argv[1])
    else:
        data = SteppingData("data/parameterSearch/kb100000000.0_kub100_expbc-0.5_t5.0_seed1.0.txt")
    print(data.rawData)
    print(data.initial_displacements)
    print(data.final_displacements)

    
