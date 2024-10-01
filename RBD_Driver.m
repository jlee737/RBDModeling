%% Ribosome Display Model

% Overview:
% 1. Bead Conjugation
% 2. IVTT Generation of stalled ribosomes
% 3. Negative Selection
% 4. Positive Selection
% 5. RNA Isolation
% 6. Reverse Transcription
% 7. PCR Amplification

%% Variables

beads.targets_pBead = 15000;
beads.pVol = 3.4 * 10^10; % Beads per mL [4]; Concentrated 4x?
beads.vol = 12.5 * 10^-6; % mL, Volume of beads used

Targets = beads.targets_pBead * beads.pVol * beads.vol;

ribo_conc = 2 * 10^-6; % M [3]

% Parameters for binding density function
pos_avg_Kd = 10^-5;
pos_std_Kd = pos_avg_Kd;

neg_avg_Kd = 10^-3;
neg_std_Kd = neg_avg_Kd;

% Washing Variables
wash_vol = 200 * 10^-6; % L
num_wash = 4;

%% IVTT Generation of stalled ribosomes

% Per 25uL reaction
IVTT.prot_conc = 8.61 * 10^-15 / (25 * 10^-6); % Roughly estimated from [4], 8.61fmol / 25uL, in M
IVTT.eff = IVTT.prot_conc / ribo_conc;

Tot_NB = 10^11;

% IVTT Complex Formation
Library.diversity = 10^11; % These variables are currently unused as diversity seems to be larger than protein amount
Library.plexity = 100; % Check with Dan

%% Negative Selection

% Equilibrium Biopanning to remove non-specific binding NBs
[Bound_NBs, tot_bound, Kd] = NegSelection(Targets * 2, Tot_NB, neg_avg_Kd, neg_std_Kd);

% Washing binders off negative selection beads
[~, tot_bound] = Washing(Tot_NB, tot_bound, Kd, Bound_NBs, wash_vol, num_wash);
Unbound_NBs = Tot_NB - tot_bound;

%% Positive Selection

% Equilbrium Biopanning to bind NBs to Targets
[Bound_NBs, tot_bound] = PosSelection(Targets, Unbound_NBs, pos_avg_Kd, pos_std_Kd);

% Washing binders off positive selection beads
[Bound_NBs, tot_bound] = Washing(Tot_NB, tot_bound, Kd, Bound_NBs, wash_vol, num_wash);

tot_bound;

%% RNA Isolation


%% Reverse Transcription


%% PCR Amplification
enrichment_fraction = tot_bound / Tot_NB;
post_PCR_dist = Bound_NBs ./ enrichment_fraction;

plot(Kd(1:100), post_PCR_dist(1:100))
title('Distribution after round of selection')

%% References

% 1. https://www.sciencedirect.com/science/article/pii/S0022283697915552#BIB22

% 2. https://pdf.sciencedirectassets.com/272314/1-s2.0-S0022519300X01981/1-s2.0-S0022519385702184/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEPb%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJGMEQCIArtGT%2FL%2BcbSB15o9Ru4YuIgBv1rnOufEg2yL8gZ2eb0AiATNwFhWIVjPwXphBLG3D8pKfQQckmbeLYuOi8lI6d7Biq8BQjv%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAUaDDA1OTAwMzU0Njg2NSIMmtEwvDe7dLcYAEvIKpAF6u41VswTxb5j0id3twmE%2FKB1Rgx%2FutssCU8DLDSZOyBc6RGNxT3JqNt73DsCq31xGkPDVHp4DDLa%2FAS4GvZByITqAapMhsTkORytzczIrMB%2Fke96SBS6iiP48WpH%2FovM7t6%2BNDKy5q6q6P6myNyV7eN2n3CLay90vjjUgcwRmahDvqADZpqzddQf6V02anYg3f61%2BTvtt0Oa5pUOtQswYxcsioR%2BdK2Rr5fDs%2FXjc15bIrPAweP46IXcyuMmpfCI6TxX05NflK5HNO2kfMmsPPHXhdU64ycRvRm4qEtNiHngdESN5cfD%2BaU5o%2FJw6VHo5k0YQ6D%2BmLA1LM8su8taQNXkZ6onP350%2BgBZ87FtdkQoXP9moOZ4hGVJL1EHOeQsrCvKMnxfA0uLQj%2BSNRFX0EjOgJc548EjlOspkX8hlzEfZziD1SKIMtQDXvfBOLZKNwSMfrCQC8DI6dSBn9tGmpLu19CxMbHJ0Lp41P1N90rvaWCHlj%2BAnCuibSAVulllpziQ1SQATw%2BKfExI4hmE65OXHJXBrPQ7jF0rwUYb3KhNo6zwm8JSw5Rd4ory6SVY3jcHSCGElQbLjLXQ6nKSS4d7V7dJQ7%2BC4YIiXRyKj%2Bs2lgXpP6zjQwE6tgv5iG3Ggk0r5zv8ZPqKHNXLF69qz3s91XnPcsr%2BiUDL%2Bkn7GapSvwT1EKENcnRchfY51n%2Bzaftrf9X1x7ylCIRoYqlMtSfHrUEhLGEOQFC9wa7I%2BWZTcDEp2Cj80uhjRhP9KNFTnKgKBVpdWrS6hxWMmJ4qclfgizWsqTHD%2F5b5anncLy1%2FC7eP0I3otCU77%2FYHYKMULovMAUSfk%2FgDRuQhklhMdixdIB9w8wLEFhtJxZZEbdkwwrX9tQY6sgE2BhZBMS%2BFWaQvlss0QhutP552iU6rrpUMnY2hvwMXRD9KA00P5Rbs6slq6d8nEUBQ39%2B89%2Bw7YVQ5tH%2F%2FZXxsYy7q0CcaCbe%2FUPuU4stwvFMGHK%2Fp0Aa2pTNgcaKMFOIHoNJXxKBIyEUFtyiiWZkWu2VdpRcC%2FMb1GJi3tCBB8VvtWepkyyOcZX3l5O7APlZ44RemL3ahb%2BPF%2FfIaiGF0wMwVLVF6qOFfzot6hoBcfvs8&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240816T142459Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY6RUDPHMW%2F20240816%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=a14cdbceb0bf64330b33a99e4a77563f1c42fc990dff756eb6de7b889a2f7d2e&hash=624eaf981440af6cf31876f4d63a52fcccf429ef6d2e38082734c7d32f8a2978&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0022519385702184&tid=spdf-8a5b59c0-332d-4c44-929b-da5f2c756948&sid=d7a28b839bbec4479f4984c4d44fd63d0077gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=15165b065657515157&rr=8b421751c8068fcc&cc=us

% 3. https://www.neb.com/en/-/media/nebus/files/pdf-faq/purexpress-faqs.pdf?rev=5e75bc74e5aa42a191a69e2b648e2b0e&hash=37662D362B00589C475256C068DFECB2

% 4. https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0015761_DynabeadsMyOneStreptavidin_T1_UG.pdf

