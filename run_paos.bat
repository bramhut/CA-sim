set /p arguments=<default_args.txt
set OMP_WAIT_POLICY=PASSIVE
set OMP_NUM_THREADS=16
out\build\x64-Release\sim-paos.exe %arguments%
PAUSE