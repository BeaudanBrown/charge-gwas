{
  inputs = {
    nixpkgs.url = "github:cachix/devenv-nixpkgs/rolling";
    systems.url = "github:nix-systems/default";
    devenv.url = "github:cachix/devenv";
    devenv.inputs.nixpkgs.follows = "nixpkgs";
  };

  nixConfig = {
    extra-trusted-public-keys = "devenv.cachix.org-1:w1cLUi8dv3hnoSPGAuibQv+f9TZLr6cv/Hm9XgU50cw=";
    extra-substituters = "https://devenv.cachix.org";
  };

  outputs = { self, nixpkgs, devenv, systems, ... } @ inputs:
    let
      forEachSystem = nixpkgs.lib.genAttrs (import systems);
    in
    {
      packages = forEachSystem (system: {
        devenv-up = self.devShells.${system}.default.config.procfileScript;
      });

      devShells = forEachSystem
        (system:
          let
            pkgs = nixpkgs.legacyPackages.${system};
          in
          {
            default = devenv.lib.mkShell {
              inherit inputs pkgs;
              modules = [
                {
                  packages = with pkgs; with rPackages; [
                    plink-ng

                    languageserver
                    glibcLocales
                    qqman
                    data_table
                    ggplot2
                    cowplot
                    dotenv
                  ];

                  languages.r = {
                    enable = true;
                  };


                  languages.python = {
                    enable = true;
                    poetry = {
                      enable = true;
                      install = {
                        enable = true;
                        installRootPackage = false;
                        onlyInstallRootPackage = false;
                        compile = false;
                        quiet = false;
                        groups = [ ];
                        ignoredGroups = [ ];
                        onlyGroups = [ ];
                        extras = [ ];
                        allExtras = false;
                        verbosity = "no";
                      };
                      activate.enable = true;
                      package = pkgs.poetry;
                    };
                  };

                  env.R_LIBS_USER = "./.Rlib";
                }
              ];
            };
          });
    };
}
