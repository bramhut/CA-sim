set /p arguments=<default_args.txt
set OMP_WAIT_POLICY=PASSIVE
out\build\x64-Release\sim-paos.exe %arguments%
PAUSE