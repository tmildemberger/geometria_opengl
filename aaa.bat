@ECHO Start of Loop

@FOR /L %%i IN (1,1,15) DO @(
  ECHO "---------" >> loglog.txt
  geometria_opengl.exe >> loglog.txt 2>&1
)