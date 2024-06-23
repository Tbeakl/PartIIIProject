import CaptureProfiling
import CaptureAttackTraces as CaptureAttack

# First do several runs to just warm up the device
for i in range(15): # 10
    print(f"Warm up {i}")
    CaptureProfiling.warm_up()

# Now actually make the profiling attacks
for i in range(85): # 85
    print(f"Profiling {i}")
    CaptureProfiling.main(i)

# Now make the attack traces where they have the same exact input and output as
# each other
for i in range(12): # 12
    print(f"Attack Counter intitial random {i}")
    CaptureAttack.main_varying_counter_initial_random(i)

for i in range(12): # 12
   print(f"Attack Counter constant {i}")
   CaptureAttack.main_constant_counter_initial_random(i)

#for i in range(11): # 12
#    print(f"Attack Counter intitial 1 {i}")
#    CaptureAttack.main_varying_counter_initial_1(i)
