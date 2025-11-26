{
  # ─────────────────────────────────── inputs ────────────────────────────────────
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    systems.url = "github:nix-systems/default";
    flake-utils = {
      url = "github:numtide/flake-utils";
      inputs.systems.follows = "systems";
    };
    pre-commit-hooks.url = "github:cachix/git-hooks.nix";
  };

  # ────────────────────────────────── outputs ───────────────────────────────────
  outputs = { self, nixpkgs, flake-utils, pre-commit-hooks, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
        {
        # ───────────────────────────────── checks ─────────────────────────────────
        checks = {
          pre-commit-check = pre-commit-hooks.lib.${system}.run {
            src = ./.;
            hooks = {
              air-fmt = {
                enable = true;
                entry = "air format";
                files = ".*\\.[rR]$";
              };
            };
          };
        };

        # ───────────────────────────────── devShells ────────────────────────────────
        devShells.default = pkgs.mkShell {
          inherit (self.checks.${system}.pre-commit-check) shellHook;
          env.R_LIBS_USER = "./.Rlib";

          buildInputs = [
            pkgs.bashInteractive
            self.checks.${system}.pre-commit-check.enabledPackages
          ];

          packages =
            with pkgs;
            [
              parallel
              plink-ng
              dos2unix
              bcftools
              samtools
              htslib
              R
              quarto
              air-formatter
              (python3.withPackages (
                python-pkgs: with python-pkgs; [
                  pandas
                  python-dotenv
                  requests
                  (buildPythonPackage rec {
                    pname = "pycap";
                    version = "2.6.0";
                    src = fetchPypi {
                      inherit pname version;
                      sha256 = "sha256-aNdAO/VzsDriTLJS+x5fc/42W2ydVMRhmQFO2v/Mj5Q=";
                    };
                    doCheck = false;
                    checkInputs = [];
                    propagatedBuildInputs = [];
                    dependencies = [semantic-version];
                  })
                ]
              ))
            ]
            ++ (with rPackages; [
              targets
              tarchetypes
              languageserver
              qqman
              data_table
              ggplot2
              cowplot
              dotenv
              tidyverse
            ]);
        };
      }
    );
}
