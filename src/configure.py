import glob

for f in glob.glob('*.cpp'):
    print('| g++ -c {} -std=c++98 -g -Werror -O2 -Wall'.format(f))

# Below are the executables

def linkme(name, objects):
    print('| g++ -o ../{} {}'.format(name, ' '.join(objects)))
    for o in objects:
        print('<', o)

linkme('generate_stepping_data', ['generate_stepping_data.o',
                                  'dynein_simulate.o',
                                  'dynein_struct_onebound.o',
                                  'dynein_struct_bothbound.o',
                                  'utilities.o',
])
linkme('simulate_unbinding_rates', ['simulate_unbinding_rates.o',
                                    'dynein_simulate.o',
                                    'dynein_struct_onebound.o',
                                    'dynein_struct_bothbound.o',
                                    'utilities.o',
])
