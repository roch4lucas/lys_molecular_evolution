#!/bin/bash

# Rodar cada modelo em paralelo
(cd M0 && codeml lysc_sites.ctl) &
(cd M1a && codeml lysc_sites.ctl) &
(cd M2a && codeml lysc_sites.ctl) &
(cd M7 && codeml lysc_sites.ctl) &
(cd M8 && codeml lysc_sites.ctl) &

# Esperar todos terminarem
wait

echo "Análises concluídas!"

