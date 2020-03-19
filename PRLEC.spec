# -*- mode: python -*-
import sys
sys.setrecursionlimit(50000)
block_cipher = None

a = Analysis(['PRLEC.py'],
             pathex=['C:\\Users\\mit\\eclipse\\workspace\\Code\\OCT_GUI'],
             binaries=[],
             datas=[],
             hiddenimports=['pywt._extensions._cwt','cython', 'sklearn', 'sklearn.neighbors.typedefs','igraph.vendor.texttable'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='PRLEC',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
	  icon='icons\\logo.ico',
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='PRLEC')
