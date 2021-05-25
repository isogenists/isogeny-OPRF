# Sage implementation of the attack against the Oblivious PRF from supersingular isogenies

The code contained in this repository can be used to win the Attack Game 12 (Auxiliary One-More SIDH) contained in the Asiacrypt 2020 [Oblivious Pseudorandom Functions from Isogenies](https://eprint.iacr.org/2020/1532) paper.

## Requirements

The code runs on Sage 9.2 compiled with Python 3.

## Usage

1. `git clone https://github.com/isogenists/isogeny-OPRF.git`
1. `cd isogeny-OPRF/src/attack`
1. `sage one_more_attack.sage`

## Credits

The `explore_dfs_optimized` has been copied from [Enric Florit](https://twitter.com/enricflorit)'s sage implementation of the [isogenies claw attack](https://gitlab.com/efz1005/isogenies)
