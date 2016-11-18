#! /usr/bin/env python2.7

import numpy, time, signal, sys, os, matplotlib
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

os.system("rm -rf PNGs") # ensure the PNGs directory is empty.
os.system("mkdir -p PNGs") # ensure the PNGs directory exists.

pe_coloring = 'energies' in sys.argv
force_vectors = 'forces' in sys.argv

view_width = 15

def close_windows(*_):
  plt.close()
  sys.exit()

usage = '''
Usage: python2 TITLE %s speed=N [show] [forces] [energies]"
       show: show animation in a window while generating movie
             omitting show makes %s faster but less exciting to watch
     forces: plot forces in movie
   energies: color circles using potential energies
''' % (sys.argv[0], sys.argv[0])

if len(sys.argv) < 2:
  print usage
  sys.exit(1)

title = sys.argv[1]

X = [0, 1, 2, 3, 4]
Y = [0, 1, 2, 3, 4]

if sys.argv[2][:6] != 'speed=':
  print usage
  print 'ERROR: second argument must start with speed=!'
  exit(1)
speed =  float(sys.argv[2][6:])

config = numpy.loadtxt("data/stepping_movie_config_" + title + ".txt")
data = numpy.genfromtxt("data/stepping_movie_data_" + title + ".txt", delimiter="\t", invalid_raise=False)
plt.ion()

if str(type(data[0])) == "<type 'numpy.float64'>":
       print "Very short animation!"
       close_windows()

gb = float(config[0]) # FIXME: make point radii based on these
gm = float(config[1])
gt = float(config[2])

ax = plt.gca()
ax.set_aspect("equal", adjustable="box")
ax.set_xlim(-view_width, view_width)
ax.set_ylim(-view_width, view_width)

microtubule = plt.plot([-view_width, view_width], [-2, -2])
plt.xlabel('$x$ (nm)')
plt.setp(microtubule, color='c', alpha=0.8, linewidth=17.0)

stalk1, = plt.plot([ X[0], X[1] ], [ Y[0], Y[1] ], color="black")
tail1,  = plt.plot([ X[1], X[2] ], [ Y[1], Y[2] ], color="black")
tail2,  = plt.plot([ X[2], X[3] ], [ Y[2], Y[3] ], color="black")
stalk2, = plt.plot([ X[3], X[4] ], [ Y[3], Y[4] ], color="black")

binding1, = plt.plot([X[0]], [Y[0]], marker='o', color="white", markersize=1)
motor1,   = plt.plot([X[1]], [Y[1]], marker='o', color="white", markersize=18)
tail,     = plt.plot([X[2]], [Y[2]], marker='o', color="red",   markersize=12)
motor2,   = plt.plot([X[3]], [Y[3]], marker='o', color="white", markersize=18)
binding2, = plt.plot([X[4]], [Y[4]], marker='o', color="white", markersize=1)

if force_vectors:
  force_line = [i for i in range(5)]
  ave_force = numpy.mean(numpy.abs(data[:,17:28]))
  force_scaling = 5 / ave_force
  for i in range(5):
    force_line[i], = plt.plot([X[i], X[i]], [Y[i],Y[i]], 'r-')

title_text = plt.text(-.9*view_width, 0.9*view_width, 'State:')
num_steps_text = plt.text(-.9*view_width, 0.8*view_width, '0 steps')
pe_text = plt.text(-view_width, 50, 'PE: ')
t_text = plt.text(-view_width+1, -view_width+1, 'Time:')

i = 0
savefigframe = 0

signal.signal(signal.SIGINT, close_windows)

num_steps = 0
previous_state = data[0][0]
while i < len(data):
  if i >= len(data):
    i = 0

  nba_scaling = min(max(1 - (data[i][2]) / (0.5*config[5]), 0), 1) # scale based on PE/kbT ratio forced between 0-1
  nma_scaling = min(max(1 - (data[i][3]) / (0.5*config[5]), 0), 1)
  ta_scaling  = min(max(1 - (data[i][4]) / (0.5*config[5]), 0), 1)
  fma_scaling = min(max(1 - (data[i][5]) / (0.5*config[5]), 0), 1)
  fba_scaling = min(max(1 - (data[i][6]) / (0.5*config[5]), 0), 1)

  if (pe_coloring == True):
    binding1.set_color((nba_scaling, 1, nba_scaling))
    motor1.set_color((nma_scaling, 1, nma_scaling))
    tail.set_color((ta_scaling, 1, ta_scaling))
    motor2.set_color((fma_scaling, 1, fma_scaling))
    binding2.set_color((fba_scaling, 1, fba_scaling))

  X[:] = data[i][7:16:2]
  Y[:] = data[i][8:17:2]
  Fx = data[i][17:27:2]
  Fy = data[i][18:28:2]

  stalk1.set_data([ X[0], X[1] ], [ Y[0], Y[1] ])
  tail1.set_data([ X[1], X[2] ], [ Y[1], Y[2] ])
  tail2.set_data([ X[2], X[3] ], [ Y[2], Y[3] ])
  stalk2.set_data([ X[3], X[4] ], [ Y[3], Y[4] ])

  if (force_vectors == True):
    for j in range(5):
      force_line[j].set_data([X[j], X[j]+force_scaling*Fx[j]], [Y[j], Y[j]+force_scaling*Fy[j]])

  binding1.set_data(X[0], Y[0])
  motor1.set_data(X[1], Y[1])
  tail.set_data(X[2], Y[2])
  motor2.set_data(X[3], Y[3])
  binding2.set_data(X[4], Y[4])

  if data[i][0] != previous_state:
    num_steps += 0.5 # we took another half step!
    previous_state = data[i][0]
    num_steps_text.set_text('%g steps' % num_steps)
  if (data[i][0] == 0):
    title_text.set_text('State: Nearbound')
    stalk1.set_linestyle('-')
    tail1.set_linestyle('-')
    tail2.set_linestyle('--')
    stalk2.set_linestyle('--')
    
    binding1.set_zorder(19)
    motor1.set_zorder(17)
    tail.set_zorder(15)
    motor2.set_zorder(13)
    binding2.set_zorder(11)

    if force_vectors:
      force_line[0].set_zorder(20)
      force_line[1].set_zorder(18)
      force_line[2].set_zorder(16)
      force_line[3].set_zorder(14)
      force_line[4].set_zorder(12)

    stalk1.set_zorder(4)
    tail1.set_zorder(3)
    tail2.set_zorder(2)
    stalk2.set_zorder(1)
    
  elif (data[i][0] == 1):
    title_text.set_text('State: Farbound')
    stalk1.set_linestyle('--')
    tail1.set_linestyle('--')
    tail2.set_linestyle('-')
    stalk2.set_linestyle('-')

    binding1.set_zorder(11)
    motor1.set_zorder(13)
    tail.set_zorder(15)
    motor2.set_zorder(17)
    binding2.set_zorder(19)

    if force_vectors:
      force_line[0].set_zorder(12)
      force_line[1].set_zorder(14)
      force_line[2].set_zorder(16)
      force_line[3].set_zorder(18)
      force_line[4].set_zorder(20)

    stalk1.set_zorder(1)
    tail1.set_zorder(2)
    tail2.set_zorder(3)
    stalk2.set_zorder(4)
    
  elif (data[i][0] == 2):
    title_text.set_text('State: Bothbound')
    stalk1.set_linestyle('-')
    tail1.set_linestyle('-')
    tail2.set_linestyle('--')
    stalk2.set_linestyle('--')

    binding1.set_zorder(19)
    motor1.set_zorder(17)
    tail.set_zorder(15)
    motor2.set_zorder(13)
    binding2.set_zorder(11)

    if force_vectors:
      force_line[0].set_zorder(20)
      force_line[1].set_zorder(18)
      force_line[2].set_zorder(16)
      force_line[3].set_zorder(14)
      force_line[4].set_zorder(12)

    stalk1.set_zorder(4)
    tail1.set_zorder(3)
    tail2.set_zorder(2)
    stalk2.set_zorder(1)
    
  elif (data[i][0] == 3):
    title_text.set_text('State: Unbound')

  pe_text.set_text('PE: %.2f' % (data[i][2]+data[i][3]+data[i][4]+data[i][5]+data[i][6]))

  # t_text.set_text("Progress: {:3.1f}%".format(i/len(data)*100))
  t_text.set_text("Time: {:g} ns".format(1e9*data[i][1]))

  i += speed
  plt.pause(0.001)
  savefigframe += 1

  fname = 'PNGs/%s-%06d.png' % (title, savefigframe)
  plt.savefig(fname)
  sys.stdout.write("video progress: %.1f%%\r" % (i/len(data)*100))
  sys.stdout.flush()

# avconv may not be present on non-Debian-related systems, in which
# case you may be able to substitute ffmpeg.
framerate = 30
avconv = "avconv -y -r %g -i PNGs/%s-%%06d.png -b 1000k movies/%s.mp4" % (framerate, title, title)
print(avconv)
os.system(avconv) # make the movie

# using convert to create a gif can be handy for small files
# os.system("convert -delay 10 PNGs/%s-*.png movies/%s.gif" % (title, title))


# os.system("mencoder -quiet PNGs/%s-*.png -mf type=png:fps=10 -ovc lavc"
#             " -lavcopts vcodec=wmv2 -oac copy -o movies/%s.mpg" % (title, title))

