import numpy as np

EPSILON = 1e-7

def equal(f1, f2):
    return abs(f1-f2) < EPSILON

class SteppingData(object):
    def __init__(self, dataFile):
        self.dataFile = dataFile
        self.rawData = np.loadtxt(self.dataFile)
        if len(self.rawData) == 8 or len(self.rawData[:,1]) <= 5:
            print("Error. File {} has less than six steps. Exiting.".format(dataFile))
            exit(1)
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

        # self.onebound_times = self.bindTimes[1:]-self.unbindTimes[1:]
        # self.bothbound_times = self.unbindTimes[1:]-self.bindTimes[:-1]
        self.onebound_times = []
        self.bothbound_times = []

        self.initial_displacements = []
        self.final_displacements = []

        self.leading_foot_steps = 0
        self.trailing_foot_steps = 0

        self.alternating_passing = 0
        self.alternating_not_passing = 0
        self.not_alternating_passing = 0
        self.not_alternating_not_passing = 0

        assert(len(self.nbx_bind)==len(self.fbx_bind))
        for s in range(1, len(self.nbx_bind)):
            if not ((self.nbx_bind[s-1] == self.nbx_bind[s]) or (self.fbx_bind[s-1] == self.fbx_bind[s])):
                print("Error, one step had both feet move")
                print(self.nbx_bind[s-1], self.nbx_bind[s], self.fbx_bind[s-1], self.fbx_bind[s])
                exit(1)
            if(self.nbx_bind[s-1]== self.nbx_bind[s] and self.fbx_bind[s-1] == self.fbx_bind[s]):
                #print("Zero-length step")
                continue
            if not (self.nbx_bind[s-1] == self.nbx_bind[s]):
                self.initial_displacements.append(self.nbx_bind[s-1]-self.fbx_bind[s-1])
                self.final_displacements.append(self.nbx_bind[s]-self.fbx_bind[s])
            elif not (self.fbx_bind[s-1] == self.fbx_bind[s]):
                self.initial_displacements.append(self.fbx_bind[s-1]-self.nbx_bind[s-1])
                self.final_displacements.append(self.fbx_bind[s] - self.nbx_bind[s])
            self.onebound_times.append(self.bindTimes[s]-self.unbindTimes[s])
            if (self.bindTimes[s] == self.unbindTimes[s]):
                print("bind and unbind are same: ", self.bindTimes[s], self.unbindTimes[s])
                exit(1)
            self.bothbound_times.append(self.unbindTimes[s]-self.bindTimes[s-1])

        for s in range(2, len(self.nbx_bind)):
            if self.nbx_bind[s-1] < self.fbx_bind[s-1]:
                self.trailing_foot = self.nbx_bind
                self.leading_foot = self.fbx_bind
            else:
                self.trailing_foot = self.fbx_bind
                self.leading_foot = self.nbx_bind
            if not equal(self.nbx_bind[s], self.nbx_bind[s-1]) and not equal(self.fbx_bind[s], self.fbx_bind[s-1]):
                print("Error, both feet moved in a step.")
                exit(1)
            if not equal(self.trailing_foot[s], self.trailing_foot[s-1]): #must've been a leading foot step
                self.leading_foot_steps += 1
                if not equal(self.trailing_foot[s-1], self.trailing_foot[s-2]): # not alternating, the last was the same foot
                    # the leading foot moved twice in a row, that makes this "not alternating"
                    self.not_alternating_not_passing += 1
                else:
                    # It is alternating because leading foot moved, but before that the other foot moved.
                    self.alternating_not_passing += 1
            else: # must've been a trailing foot step
                self.trailing_foot_steps += 1
                if equal(self.trailing_foot[s-1], self.trailing_foot[s-2]): # alternating, other foot moved last time
                    if self.trailing_foot[s] > self.leading_foot[s]:
                        self.not_alternating_passing += 1
                    else:
                        self.not_alternating_not_passing += 1
                else:
                    if self.trailing_foot[s] > self.leading_foot[s]:
                        self.alternating_passing += 1
                    else:
                        self.alternating_not_passing += 1

        self.initial_displacements = np.asarray(self.initial_displacements)
        self.final_displacements = np.asarray(self.final_displacements)

        self.step_lengths = self.final_displacements - self.initial_displacements
        self.step_times = self.onebound_times + self.bothbound_times



class MovieData(object):
    def __init__(self, movieFile):
        self.rawData = np.loadtxt(movieFile)
        self.times = self.rawData[:, 1]
        self.PE1 = self.rawData[:, 2]
        self.PE2 = self.rawData[:, 3]
        self.PE3 = self.rawData[:, 4]
        self.PE4 = self.rawData[:, 5]
        self.PE5 = self.rawData[:, 6]
        self.x1 = self.rawData[:, 7]
        self.y1 = self.rawData[:, 8]
        self.x2 = self.rawData[:, 9]
        self.y2 = self.rawData[:, 10]
        self.x3 = self.rawData[:, 11]
        self.y3 = self.rawData[:, 12]
        self.x4 = self.rawData[:, 13]
        self.y4 = self.rawData[:, 14]
        self.x5 = self.rawData[:, 15]
        self.y5 = self.rawData[:, 16]
        self.fx1 = self.rawData[:, 17]
        self.fy1 = self.rawData[:, 18]
        self.fx2 = self.rawData[:, 19]
        self.fy2 = self.rawData[:, 20]
        self.fx3 = self.rawData[:, 21]
        self.fy3 = self.rawData[:, 22]
        self.fx4 = self.rawData[:, 23]
        self.fy4 = self.rawData[:, 24]
        self.fx5 = self.rawData[:, 25]
        self.fy5 = self.rawData[:, 26]


if __name__ == "__main__":
    data = MovieData('/home/john/gitRepos/dynein_walk/data/paper_trajectory_movie_data.txt')
    print(np.shape(data.rawData))
    print(np.shape(data.times))




























