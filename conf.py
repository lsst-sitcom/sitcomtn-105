# See the Documenteer docs for how to customize conf.py:
# https://documenteer.lsst.io/technotes/

from documenteer.conf.technote import *  # noqa F401 F403

extensions += ['nbsphinx']
nbsphinx_execute = 'never'
exclude_patterns.append("notebooks")
