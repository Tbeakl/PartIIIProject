import chipwhisperer as cw
import time 
scope = cw.scope()

target = cw.target(scope, cw.targets.SimpleSerial) #cw.targets.SimpleSerial can be omitted
scope.default_setup()
cw.program_target(scope, cw.programmers.STM32FProgrammer, "simpleserial-base-CWLITEARM-volatile.hex")
#cw.program_target(scope, cw.programmers.STM32FProgrammer, "simpleserial-8-on-32-CWLITEARM.hex")
