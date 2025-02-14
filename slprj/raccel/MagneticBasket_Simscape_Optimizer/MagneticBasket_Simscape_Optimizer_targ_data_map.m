    function targMap = targDataMap(),

    ;%***********************
    ;% Create Parameter Map *
    ;%***********************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 3;
        sectIdxOffset = 0;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc paramMap
        ;%
        paramMap.nSections           = nTotSects;
        paramMap.sectIdxOffset       = sectIdxOffset;
            paramMap.sections(nTotSects) = dumSection; %prealloc
        paramMap.nTotData            = -1;

        ;%
        ;% Auto data (rtP)
        ;%
            section.nData     = 1;
            section.data(1)  = dumData; %prealloc

                    ;% rtP.Fixed
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

            nTotData = nTotData + section.nData;
            paramMap.sections(1) = section;
            clear section

            section.nData     = 1;
            section.data(1)  = dumData; %prealloc

                    ;% rtP.SimParams
                    section.data(1).logicalSrcIdx = 1;
                    section.data(1).dtTransOffset = 0;

            nTotData = nTotData + section.nData;
            paramMap.sections(2) = section;
            clear section

            section.nData     = 101;
            section.data(101)  = dumData; %prealloc

                    ;% rtP.x
                    section.data(1).logicalSrcIdx = 2;
                    section.data(1).dtTransOffset = 0;

                    ;% rtP.NormalizeVector_maxzero
                    section.data(2).logicalSrcIdx = 3;
                    section.data(2).dtTransOffset = 14;

                    ;% rtP.NormalizeVector_maxzero_ab1bjfgz5c
                    section.data(3).logicalSrcIdx = 4;
                    section.data(3).dtTransOffset = 15;

                    ;% rtP.NormalizeVector_maxzero_jpxndxsk21
                    section.data(4).logicalSrcIdx = 5;
                    section.data(4).dtTransOffset = 16;

                    ;% rtP.NormalizeVector1_maxzero
                    section.data(5).logicalSrcIdx = 6;
                    section.data(5).dtTransOffset = 17;

                    ;% rtP.NormalizeVector_maxzero_hach1vuaw5
                    section.data(6).logicalSrcIdx = 7;
                    section.data(6).dtTransOffset = 18;

                    ;% rtP.NormalizeVector_maxzero_nzh4nf0nvr
                    section.data(7).logicalSrcIdx = 8;
                    section.data(7).dtTransOffset = 19;

                    ;% rtP.NormalizeVector_maxzero_pjti1ch2qi
                    section.data(8).logicalSrcIdx = 9;
                    section.data(8).dtTransOffset = 20;

                    ;% rtP.NormalizeVector1_maxzero_hfyr12og40
                    section.data(9).logicalSrcIdx = 10;
                    section.data(9).dtTransOffset = 21;

                    ;% rtP.NormalizeVector_maxzero_jvlgwkux2b
                    section.data(10).logicalSrcIdx = 11;
                    section.data(10).dtTransOffset = 22;

                    ;% rtP.NormalizeVector_maxzero_k3zfgywmre
                    section.data(11).logicalSrcIdx = 12;
                    section.data(11).dtTransOffset = 23;

                    ;% rtP.NormalizeVector_maxzero_jceokergvn
                    section.data(12).logicalSrcIdx = 13;
                    section.data(12).dtTransOffset = 24;

                    ;% rtP.NormalizeVector1_maxzero_edxcmuydlr
                    section.data(13).logicalSrcIdx = 14;
                    section.data(13).dtTransOffset = 25;

                    ;% rtP.NormalizeVector_maxzero_levzkkbbsj
                    section.data(14).logicalSrcIdx = 15;
                    section.data(14).dtTransOffset = 26;

                    ;% rtP.NormalizeVector_maxzero_kt51ehrppa
                    section.data(15).logicalSrcIdx = 16;
                    section.data(15).dtTransOffset = 27;

                    ;% rtP.NormalizeVector_maxzero_hbfzo5qsh2
                    section.data(16).logicalSrcIdx = 17;
                    section.data(16).dtTransOffset = 28;

                    ;% rtP.NormalizeVector1_maxzero_os1qknjcxt
                    section.data(17).logicalSrcIdx = 18;
                    section.data(17).dtTransOffset = 29;

                    ;% rtP.NormalizeVector_maxzero_anlse3xi4h
                    section.data(18).logicalSrcIdx = 19;
                    section.data(18).dtTransOffset = 30;

                    ;% rtP.NormalizeVector_maxzero_dfin15pekf
                    section.data(19).logicalSrcIdx = 20;
                    section.data(19).dtTransOffset = 31;

                    ;% rtP.NormalizeVector_maxzero_nfiirapz00
                    section.data(20).logicalSrcIdx = 21;
                    section.data(20).dtTransOffset = 32;

                    ;% rtP.NormalizeVector1_maxzero_pbowfqqdx5
                    section.data(21).logicalSrcIdx = 22;
                    section.data(21).dtTransOffset = 33;

                    ;% rtP.NormalizeVector_maxzero_ap3h53mgxx
                    section.data(22).logicalSrcIdx = 23;
                    section.data(22).dtTransOffset = 34;

                    ;% rtP.NormalizeVector_maxzero_kytveykcqy
                    section.data(23).logicalSrcIdx = 24;
                    section.data(23).dtTransOffset = 35;

                    ;% rtP.NormalizeVector_maxzero_ihmtmp1vz1
                    section.data(24).logicalSrcIdx = 25;
                    section.data(24).dtTransOffset = 36;

                    ;% rtP.NormalizeVector1_maxzero_byir5pwee4
                    section.data(25).logicalSrcIdx = 26;
                    section.data(25).dtTransOffset = 37;

                    ;% rtP.NormalizeVector_maxzero_fqwaltgljo
                    section.data(26).logicalSrcIdx = 27;
                    section.data(26).dtTransOffset = 38;

                    ;% rtP.NormalizeVector_maxzero_hkbi1qbiut
                    section.data(27).logicalSrcIdx = 28;
                    section.data(27).dtTransOffset = 39;

                    ;% rtP.NormalizeVector_maxzero_dnlmszn2ri
                    section.data(28).logicalSrcIdx = 29;
                    section.data(28).dtTransOffset = 40;

                    ;% rtP.NormalizeVector1_maxzero_j5bp1hk5bd
                    section.data(29).logicalSrcIdx = 30;
                    section.data(29).dtTransOffset = 41;

                    ;% rtP.Gain_Gain
                    section.data(30).logicalSrcIdx = 31;
                    section.data(30).dtTransOffset = 42;

                    ;% rtP.Gain_Gain_iby0aikld1
                    section.data(31).logicalSrcIdx = 32;
                    section.data(31).dtTransOffset = 43;

                    ;% rtP.Gain1_Gain
                    section.data(32).logicalSrcIdx = 33;
                    section.data(32).dtTransOffset = 44;

                    ;% rtP.Gain_Gain_o3ussx2pey
                    section.data(33).logicalSrcIdx = 34;
                    section.data(33).dtTransOffset = 45;

                    ;% rtP.Gain1_Gain_inzuw1ftmc
                    section.data(34).logicalSrcIdx = 35;
                    section.data(34).dtTransOffset = 46;

                    ;% rtP.Gain_Gain_ch2juginai
                    section.data(35).logicalSrcIdx = 36;
                    section.data(35).dtTransOffset = 47;

                    ;% rtP.Gain1_Gain_iijf52nd0p
                    section.data(36).logicalSrcIdx = 37;
                    section.data(36).dtTransOffset = 48;

                    ;% rtP.Gain_Gain_cnz04octod
                    section.data(37).logicalSrcIdx = 38;
                    section.data(37).dtTransOffset = 49;

                    ;% rtP.Gain1_Gain_kgqmwgbxmo
                    section.data(38).logicalSrcIdx = 39;
                    section.data(38).dtTransOffset = 50;

                    ;% rtP.Gain_Gain_mxfyehzrbo
                    section.data(39).logicalSrcIdx = 40;
                    section.data(39).dtTransOffset = 51;

                    ;% rtP.Gain1_Gain_pklggewqk2
                    section.data(40).logicalSrcIdx = 41;
                    section.data(40).dtTransOffset = 52;

                    ;% rtP.Gain_Gain_ch2qj2nra0
                    section.data(41).logicalSrcIdx = 42;
                    section.data(41).dtTransOffset = 53;

                    ;% rtP.Gain1_Gain_mt5ivn2ky5
                    section.data(42).logicalSrcIdx = 43;
                    section.data(42).dtTransOffset = 54;

                    ;% rtP.Gain_Gain_kzgbyxe2u1
                    section.data(43).logicalSrcIdx = 44;
                    section.data(43).dtTransOffset = 55;

                    ;% rtP.Gain1_Gain_fzh3lvxzcm
                    section.data(44).logicalSrcIdx = 45;
                    section.data(44).dtTransOffset = 56;

                    ;% rtP.Gain_Gain_gnnygtgrhf
                    section.data(45).logicalSrcIdx = 46;
                    section.data(45).dtTransOffset = 57;

                    ;% rtP.Constant_Value
                    section.data(46).logicalSrcIdx = 47;
                    section.data(46).dtTransOffset = 58;

                    ;% rtP.Constant_Value_foidw3bu40
                    section.data(47).logicalSrcIdx = 48;
                    section.data(47).dtTransOffset = 59;

                    ;% rtP.MagnetDipoleMoment_Value
                    section.data(48).logicalSrcIdx = 49;
                    section.data(48).dtTransOffset = 60;

                    ;% rtP.Constant_Value_ffpsqzlbdy
                    section.data(49).logicalSrcIdx = 50;
                    section.data(49).dtTransOffset = 63;

                    ;% rtP.Constant1_Value
                    section.data(50).logicalSrcIdx = 51;
                    section.data(50).dtTransOffset = 64;

                    ;% rtP.Constant_Value_aldt0rz0ce
                    section.data(51).logicalSrcIdx = 52;
                    section.data(51).dtTransOffset = 73;

                    ;% rtP.Constant_Value_hsgo25r4sk
                    section.data(52).logicalSrcIdx = 53;
                    section.data(52).dtTransOffset = 74;

                    ;% rtP.Constant_Value_lu1rr1djht
                    section.data(53).logicalSrcIdx = 54;
                    section.data(53).dtTransOffset = 75;

                    ;% rtP.Constant_Value_hrhqyti5xa
                    section.data(54).logicalSrcIdx = 55;
                    section.data(54).dtTransOffset = 76;

                    ;% rtP.Constant_Value_bgbfcrykss
                    section.data(55).logicalSrcIdx = 56;
                    section.data(55).dtTransOffset = 77;

                    ;% rtP.MagnetDipoleMoment_Value_hus3icj2zx
                    section.data(56).logicalSrcIdx = 57;
                    section.data(56).dtTransOffset = 78;

                    ;% rtP.Constant_Value_ccdbgnkv2m
                    section.data(57).logicalSrcIdx = 58;
                    section.data(57).dtTransOffset = 81;

                    ;% rtP.Constant1_Value_mat3dcwilq
                    section.data(58).logicalSrcIdx = 59;
                    section.data(58).dtTransOffset = 82;

                    ;% rtP.Constant_Value_ag31g0q0wg
                    section.data(59).logicalSrcIdx = 60;
                    section.data(59).dtTransOffset = 91;

                    ;% rtP.Constant_Value_lpo3xrcqjh
                    section.data(60).logicalSrcIdx = 61;
                    section.data(60).dtTransOffset = 92;

                    ;% rtP.Constant_Value_kebaex1duz
                    section.data(61).logicalSrcIdx = 62;
                    section.data(61).dtTransOffset = 93;

                    ;% rtP.Constant_Value_crcwx2inuv
                    section.data(62).logicalSrcIdx = 63;
                    section.data(62).dtTransOffset = 94;

                    ;% rtP.Constant_Value_p4qexwk3wd
                    section.data(63).logicalSrcIdx = 64;
                    section.data(63).dtTransOffset = 95;

                    ;% rtP.MagnetDipoleMoment_Value_m4h5ot0d2d
                    section.data(64).logicalSrcIdx = 65;
                    section.data(64).dtTransOffset = 96;

                    ;% rtP.Constant_Value_pxe2jyk15f
                    section.data(65).logicalSrcIdx = 66;
                    section.data(65).dtTransOffset = 99;

                    ;% rtP.Constant1_Value_obh4dra0ft
                    section.data(66).logicalSrcIdx = 67;
                    section.data(66).dtTransOffset = 100;

                    ;% rtP.Constant_Value_a4fqp5ywdt
                    section.data(67).logicalSrcIdx = 68;
                    section.data(67).dtTransOffset = 109;

                    ;% rtP.Constant_Value_lwhq5q2fck
                    section.data(68).logicalSrcIdx = 69;
                    section.data(68).dtTransOffset = 110;

                    ;% rtP.Constant_Value_mqmlcee2pu
                    section.data(69).logicalSrcIdx = 70;
                    section.data(69).dtTransOffset = 111;

                    ;% rtP.Constant_Value_mulpe1a5eh
                    section.data(70).logicalSrcIdx = 71;
                    section.data(70).dtTransOffset = 112;

                    ;% rtP.Constant_Value_f4bgccfy2v
                    section.data(71).logicalSrcIdx = 72;
                    section.data(71).dtTransOffset = 113;

                    ;% rtP.MagnetDipoleMoment_Value_am3fjqdljc
                    section.data(72).logicalSrcIdx = 73;
                    section.data(72).dtTransOffset = 114;

                    ;% rtP.Constant_Value_kyzvzubrmq
                    section.data(73).logicalSrcIdx = 74;
                    section.data(73).dtTransOffset = 117;

                    ;% rtP.Constant1_Value_htaifltrf1
                    section.data(74).logicalSrcIdx = 75;
                    section.data(74).dtTransOffset = 118;

                    ;% rtP.Constant_Value_hooldt5mbh
                    section.data(75).logicalSrcIdx = 76;
                    section.data(75).dtTransOffset = 127;

                    ;% rtP.Constant_Value_fyrymlyvwp
                    section.data(76).logicalSrcIdx = 77;
                    section.data(76).dtTransOffset = 128;

                    ;% rtP.Constant_Value_iinhm22lbs
                    section.data(77).logicalSrcIdx = 78;
                    section.data(77).dtTransOffset = 129;

                    ;% rtP.Constant_Value_o54qqmia1e
                    section.data(78).logicalSrcIdx = 79;
                    section.data(78).dtTransOffset = 130;

                    ;% rtP.Constant_Value_fkwr40m5ih
                    section.data(79).logicalSrcIdx = 80;
                    section.data(79).dtTransOffset = 131;

                    ;% rtP.MagnetDipoleMoment_Value_hqiqmgkprz
                    section.data(80).logicalSrcIdx = 81;
                    section.data(80).dtTransOffset = 132;

                    ;% rtP.Constant_Value_c50rewkplm
                    section.data(81).logicalSrcIdx = 82;
                    section.data(81).dtTransOffset = 135;

                    ;% rtP.Constant1_Value_ajigyq4ohn
                    section.data(82).logicalSrcIdx = 83;
                    section.data(82).dtTransOffset = 136;

                    ;% rtP.Constant_Value_gxy2qdgtgj
                    section.data(83).logicalSrcIdx = 84;
                    section.data(83).dtTransOffset = 145;

                    ;% rtP.Constant_Value_oqtflizduh
                    section.data(84).logicalSrcIdx = 85;
                    section.data(84).dtTransOffset = 146;

                    ;% rtP.Constant_Value_e34xedeirm
                    section.data(85).logicalSrcIdx = 86;
                    section.data(85).dtTransOffset = 147;

                    ;% rtP.Constant_Value_l5p2q1mni4
                    section.data(86).logicalSrcIdx = 87;
                    section.data(86).dtTransOffset = 148;

                    ;% rtP.Constant_Value_cuvlenn5ze
                    section.data(87).logicalSrcIdx = 88;
                    section.data(87).dtTransOffset = 149;

                    ;% rtP.MagnetDipoleMoment_Value_clffg1mnnq
                    section.data(88).logicalSrcIdx = 89;
                    section.data(88).dtTransOffset = 150;

                    ;% rtP.Constant_Value_pik4yifivl
                    section.data(89).logicalSrcIdx = 90;
                    section.data(89).dtTransOffset = 153;

                    ;% rtP.Constant1_Value_dagdlneepy
                    section.data(90).logicalSrcIdx = 91;
                    section.data(90).dtTransOffset = 154;

                    ;% rtP.Constant_Value_ndmgh4xmik
                    section.data(91).logicalSrcIdx = 92;
                    section.data(91).dtTransOffset = 163;

                    ;% rtP.Constant_Value_kqizw1abe3
                    section.data(92).logicalSrcIdx = 93;
                    section.data(92).dtTransOffset = 164;

                    ;% rtP.Constant_Value_gimcsalf5s
                    section.data(93).logicalSrcIdx = 94;
                    section.data(93).dtTransOffset = 165;

                    ;% rtP.Constant_Value_h2tkrf03re
                    section.data(94).logicalSrcIdx = 95;
                    section.data(94).dtTransOffset = 166;

                    ;% rtP.Constant_Value_oc3ljvgxhv
                    section.data(95).logicalSrcIdx = 96;
                    section.data(95).dtTransOffset = 167;

                    ;% rtP.MagnetDipoleMoment_Value_af4j3ouf2p
                    section.data(96).logicalSrcIdx = 97;
                    section.data(96).dtTransOffset = 168;

                    ;% rtP.Constant_Value_nu3mnt2521
                    section.data(97).logicalSrcIdx = 98;
                    section.data(97).dtTransOffset = 171;

                    ;% rtP.Constant1_Value_k2adofcs3i
                    section.data(98).logicalSrcIdx = 99;
                    section.data(98).dtTransOffset = 172;

                    ;% rtP.Constant_Value_m0tx3vya3h
                    section.data(99).logicalSrcIdx = 100;
                    section.data(99).dtTransOffset = 181;

                    ;% rtP.Constant_Value_jgubgdhktv
                    section.data(100).logicalSrcIdx = 101;
                    section.data(100).dtTransOffset = 182;

                    ;% rtP.Constant_Value_nfejvnig4t
                    section.data(101).logicalSrcIdx = 102;
                    section.data(101).dtTransOffset = 183;

            nTotData = nTotData + section.nData;
            paramMap.sections(3) = section;
            clear section


            ;%
            ;% Non-auto Data (parameter)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        paramMap.nTotData = nTotData;



    ;%**************************
    ;% Create Block Output Map *
    ;%**************************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 1;
        sectIdxOffset = 0;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc sigMap
        ;%
        sigMap.nSections           = nTotSects;
        sigMap.sectIdxOffset       = sectIdxOffset;
            sigMap.sections(nTotSects) = dumSection; %prealloc
        sigMap.nTotData            = -1;

        ;%
        ;% Auto data (rtB)
        ;%
            section.nData     = 333;
            section.data(333)  = dumData; %prealloc

                    ;% rtB.aay54cmgia
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% rtB.nfyrk0gene
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 18;

                    ;% rtB.gdleddk2bu
                    section.data(3).logicalSrcIdx = 2;
                    section.data(3).dtTransOffset = 172;

                    ;% rtB.o3cnonqf2y
                    section.data(4).logicalSrcIdx = 3;
                    section.data(4).dtTransOffset = 179;

                    ;% rtB.bjyrxtyjst
                    section.data(5).logicalSrcIdx = 4;
                    section.data(5).dtTransOffset = 186;

                    ;% rtB.ecnsyysan1
                    section.data(6).logicalSrcIdx = 5;
                    section.data(6).dtTransOffset = 189;

                    ;% rtB.d11srgkxv1
                    section.data(7).logicalSrcIdx = 6;
                    section.data(7).dtTransOffset = 192;

                    ;% rtB.oywwaak55z
                    section.data(8).logicalSrcIdx = 7;
                    section.data(8).dtTransOffset = 195;

                    ;% rtB.n3uyejqyiv
                    section.data(9).logicalSrcIdx = 8;
                    section.data(9).dtTransOffset = 196;

                    ;% rtB.dg2tz20qh3
                    section.data(10).logicalSrcIdx = 9;
                    section.data(10).dtTransOffset = 199;

                    ;% rtB.jd0epngxde
                    section.data(11).logicalSrcIdx = 10;
                    section.data(11).dtTransOffset = 208;

                    ;% rtB.i1wwep0quu
                    section.data(12).logicalSrcIdx = 11;
                    section.data(12).dtTransOffset = 209;

                    ;% rtB.jytdis5y4s
                    section.data(13).logicalSrcIdx = 12;
                    section.data(13).dtTransOffset = 210;

                    ;% rtB.dvlkqx1wrk
                    section.data(14).logicalSrcIdx = 13;
                    section.data(14).dtTransOffset = 214;

                    ;% rtB.axlew0gacz
                    section.data(15).logicalSrcIdx = 14;
                    section.data(15).dtTransOffset = 217;

                    ;% rtB.ffk0uso1jv
                    section.data(16).logicalSrcIdx = 15;
                    section.data(16).dtTransOffset = 220;

                    ;% rtB.j2caeelflf
                    section.data(17).logicalSrcIdx = 16;
                    section.data(17).dtTransOffset = 229;

                    ;% rtB.gzdxak0rgb
                    section.data(18).logicalSrcIdx = 17;
                    section.data(18).dtTransOffset = 238;

                    ;% rtB.nhpvoyt4rs
                    section.data(19).logicalSrcIdx = 18;
                    section.data(19).dtTransOffset = 247;

                    ;% rtB.mnl5wynqkw
                    section.data(20).logicalSrcIdx = 19;
                    section.data(20).dtTransOffset = 250;

                    ;% rtB.n5zmerun4k
                    section.data(21).logicalSrcIdx = 20;
                    section.data(21).dtTransOffset = 251;

                    ;% rtB.i210rcby1e
                    section.data(22).logicalSrcIdx = 21;
                    section.data(22).dtTransOffset = 254;

                    ;% rtB.njlotg3eyk
                    section.data(23).logicalSrcIdx = 22;
                    section.data(23).dtTransOffset = 260;

                    ;% rtB.ciemy5f0t3
                    section.data(24).logicalSrcIdx = 23;
                    section.data(24).dtTransOffset = 263;

                    ;% rtB.i3fsmpf0l4
                    section.data(25).logicalSrcIdx = 24;
                    section.data(25).dtTransOffset = 266;

                    ;% rtB.p1t00nt2c4
                    section.data(26).logicalSrcIdx = 25;
                    section.data(26).dtTransOffset = 269;

                    ;% rtB.lkk1p3hkce
                    section.data(27).logicalSrcIdx = 26;
                    section.data(27).dtTransOffset = 270;

                    ;% rtB.j0nzmiwajs
                    section.data(28).logicalSrcIdx = 27;
                    section.data(28).dtTransOffset = 271;

                    ;% rtB.p4b33uxmxl
                    section.data(29).logicalSrcIdx = 28;
                    section.data(29).dtTransOffset = 272;

                    ;% rtB.kerf1yzlz0
                    section.data(30).logicalSrcIdx = 29;
                    section.data(30).dtTransOffset = 276;

                    ;% rtB.dxk5dr2guu
                    section.data(31).logicalSrcIdx = 30;
                    section.data(31).dtTransOffset = 279;

                    ;% rtB.gfwkvsxkxh
                    section.data(32).logicalSrcIdx = 31;
                    section.data(32).dtTransOffset = 282;

                    ;% rtB.jmgra3bzio
                    section.data(33).logicalSrcIdx = 32;
                    section.data(33).dtTransOffset = 285;

                    ;% rtB.n43m55vtpu
                    section.data(34).logicalSrcIdx = 33;
                    section.data(34).dtTransOffset = 286;

                    ;% rtB.dm0nsvocoi
                    section.data(35).logicalSrcIdx = 34;
                    section.data(35).dtTransOffset = 287;

                    ;% rtB.pzoqxgrplf
                    section.data(36).logicalSrcIdx = 35;
                    section.data(36).dtTransOffset = 290;

                    ;% rtB.pi3p3kites
                    section.data(37).logicalSrcIdx = 36;
                    section.data(37).dtTransOffset = 293;

                    ;% rtB.j4tq1cqvap
                    section.data(38).logicalSrcIdx = 37;
                    section.data(38).dtTransOffset = 296;

                    ;% rtB.d4sl1kvckp
                    section.data(39).logicalSrcIdx = 38;
                    section.data(39).dtTransOffset = 299;

                    ;% rtB.k0e5tb5t4r
                    section.data(40).logicalSrcIdx = 39;
                    section.data(40).dtTransOffset = 300;

                    ;% rtB.mrse3gkxz4
                    section.data(41).logicalSrcIdx = 40;
                    section.data(41).dtTransOffset = 303;

                    ;% rtB.hjkncigkim
                    section.data(42).logicalSrcIdx = 41;
                    section.data(42).dtTransOffset = 306;

                    ;% rtB.bv20qr1sq1
                    section.data(43).logicalSrcIdx = 42;
                    section.data(43).dtTransOffset = 309;

                    ;% rtB.jknqpcyts5
                    section.data(44).logicalSrcIdx = 43;
                    section.data(44).dtTransOffset = 312;

                    ;% rtB.fo44k43jik
                    section.data(45).logicalSrcIdx = 44;
                    section.data(45).dtTransOffset = 313;

                    ;% rtB.leoubbnbua
                    section.data(46).logicalSrcIdx = 45;
                    section.data(46).dtTransOffset = 316;

                    ;% rtB.eb5u3a2x55
                    section.data(47).logicalSrcIdx = 46;
                    section.data(47).dtTransOffset = 325;

                    ;% rtB.eded4se0e2
                    section.data(48).logicalSrcIdx = 47;
                    section.data(48).dtTransOffset = 326;

                    ;% rtB.epqt4glrc2
                    section.data(49).logicalSrcIdx = 48;
                    section.data(49).dtTransOffset = 327;

                    ;% rtB.huyyr23tvs
                    section.data(50).logicalSrcIdx = 49;
                    section.data(50).dtTransOffset = 331;

                    ;% rtB.a5iwna5ofw
                    section.data(51).logicalSrcIdx = 50;
                    section.data(51).dtTransOffset = 334;

                    ;% rtB.n5qkyzyc0p
                    section.data(52).logicalSrcIdx = 51;
                    section.data(52).dtTransOffset = 337;

                    ;% rtB.cviiadazje
                    section.data(53).logicalSrcIdx = 52;
                    section.data(53).dtTransOffset = 346;

                    ;% rtB.hrdkh50xl2
                    section.data(54).logicalSrcIdx = 53;
                    section.data(54).dtTransOffset = 355;

                    ;% rtB.p4bi154apu
                    section.data(55).logicalSrcIdx = 54;
                    section.data(55).dtTransOffset = 364;

                    ;% rtB.ollqza00rh
                    section.data(56).logicalSrcIdx = 55;
                    section.data(56).dtTransOffset = 367;

                    ;% rtB.kypr3ufomq
                    section.data(57).logicalSrcIdx = 56;
                    section.data(57).dtTransOffset = 368;

                    ;% rtB.ljklcjruel
                    section.data(58).logicalSrcIdx = 57;
                    section.data(58).dtTransOffset = 371;

                    ;% rtB.giyavk0qxh
                    section.data(59).logicalSrcIdx = 58;
                    section.data(59).dtTransOffset = 377;

                    ;% rtB.iwp0nxhxfg
                    section.data(60).logicalSrcIdx = 59;
                    section.data(60).dtTransOffset = 380;

                    ;% rtB.cwwkrda1qf
                    section.data(61).logicalSrcIdx = 60;
                    section.data(61).dtTransOffset = 383;

                    ;% rtB.ojzjawyp2p
                    section.data(62).logicalSrcIdx = 61;
                    section.data(62).dtTransOffset = 386;

                    ;% rtB.j22ln0sz4g
                    section.data(63).logicalSrcIdx = 62;
                    section.data(63).dtTransOffset = 387;

                    ;% rtB.ltesgv5ex5
                    section.data(64).logicalSrcIdx = 63;
                    section.data(64).dtTransOffset = 388;

                    ;% rtB.e4k3bfcn1z
                    section.data(65).logicalSrcIdx = 64;
                    section.data(65).dtTransOffset = 389;

                    ;% rtB.fi20a5zfv5
                    section.data(66).logicalSrcIdx = 65;
                    section.data(66).dtTransOffset = 393;

                    ;% rtB.oexwmjcbqo
                    section.data(67).logicalSrcIdx = 66;
                    section.data(67).dtTransOffset = 396;

                    ;% rtB.ljdiomnr0e
                    section.data(68).logicalSrcIdx = 67;
                    section.data(68).dtTransOffset = 399;

                    ;% rtB.nj4zxwlnaz
                    section.data(69).logicalSrcIdx = 68;
                    section.data(69).dtTransOffset = 402;

                    ;% rtB.b4eypdlrbg
                    section.data(70).logicalSrcIdx = 69;
                    section.data(70).dtTransOffset = 403;

                    ;% rtB.mswvoafvyk
                    section.data(71).logicalSrcIdx = 70;
                    section.data(71).dtTransOffset = 404;

                    ;% rtB.ga3a4xdbbc
                    section.data(72).logicalSrcIdx = 71;
                    section.data(72).dtTransOffset = 407;

                    ;% rtB.eecvbm2tyj
                    section.data(73).logicalSrcIdx = 72;
                    section.data(73).dtTransOffset = 410;

                    ;% rtB.bkxjvl2wwl
                    section.data(74).logicalSrcIdx = 73;
                    section.data(74).dtTransOffset = 413;

                    ;% rtB.jcdegmljub
                    section.data(75).logicalSrcIdx = 74;
                    section.data(75).dtTransOffset = 414;

                    ;% rtB.dakqbtmvna
                    section.data(76).logicalSrcIdx = 75;
                    section.data(76).dtTransOffset = 417;

                    ;% rtB.j54ujtxdmj
                    section.data(77).logicalSrcIdx = 76;
                    section.data(77).dtTransOffset = 420;

                    ;% rtB.ol0mge5uht
                    section.data(78).logicalSrcIdx = 77;
                    section.data(78).dtTransOffset = 423;

                    ;% rtB.fusdzrgc4y
                    section.data(79).logicalSrcIdx = 78;
                    section.data(79).dtTransOffset = 426;

                    ;% rtB.n2cw123fsj
                    section.data(80).logicalSrcIdx = 79;
                    section.data(80).dtTransOffset = 429;

                    ;% rtB.gtc4x0vmrj
                    section.data(81).logicalSrcIdx = 80;
                    section.data(81).dtTransOffset = 430;

                    ;% rtB.igdgu1b2hk
                    section.data(82).logicalSrcIdx = 81;
                    section.data(82).dtTransOffset = 433;

                    ;% rtB.lntlwwhs4g
                    section.data(83).logicalSrcIdx = 82;
                    section.data(83).dtTransOffset = 442;

                    ;% rtB.lvldddhyjb
                    section.data(84).logicalSrcIdx = 83;
                    section.data(84).dtTransOffset = 443;

                    ;% rtB.j1izfy2apg
                    section.data(85).logicalSrcIdx = 84;
                    section.data(85).dtTransOffset = 444;

                    ;% rtB.f2e0m2twoo
                    section.data(86).logicalSrcIdx = 85;
                    section.data(86).dtTransOffset = 448;

                    ;% rtB.pctih4a2gv
                    section.data(87).logicalSrcIdx = 86;
                    section.data(87).dtTransOffset = 451;

                    ;% rtB.d15mjte2mm
                    section.data(88).logicalSrcIdx = 87;
                    section.data(88).dtTransOffset = 454;

                    ;% rtB.o0sl3javac
                    section.data(89).logicalSrcIdx = 88;
                    section.data(89).dtTransOffset = 463;

                    ;% rtB.lk2g50sct1
                    section.data(90).logicalSrcIdx = 89;
                    section.data(90).dtTransOffset = 472;

                    ;% rtB.kucavif4r4
                    section.data(91).logicalSrcIdx = 90;
                    section.data(91).dtTransOffset = 481;

                    ;% rtB.fd02dku0jj
                    section.data(92).logicalSrcIdx = 91;
                    section.data(92).dtTransOffset = 484;

                    ;% rtB.j3daqzj5si
                    section.data(93).logicalSrcIdx = 92;
                    section.data(93).dtTransOffset = 485;

                    ;% rtB.e30x0ooogi
                    section.data(94).logicalSrcIdx = 93;
                    section.data(94).dtTransOffset = 488;

                    ;% rtB.dgv155wtql
                    section.data(95).logicalSrcIdx = 94;
                    section.data(95).dtTransOffset = 494;

                    ;% rtB.bcueeadxdi
                    section.data(96).logicalSrcIdx = 95;
                    section.data(96).dtTransOffset = 497;

                    ;% rtB.nalsmrc4um
                    section.data(97).logicalSrcIdx = 96;
                    section.data(97).dtTransOffset = 500;

                    ;% rtB.gmj04jp0xd
                    section.data(98).logicalSrcIdx = 97;
                    section.data(98).dtTransOffset = 503;

                    ;% rtB.ij2rjc05nz
                    section.data(99).logicalSrcIdx = 98;
                    section.data(99).dtTransOffset = 504;

                    ;% rtB.pyxwiqcuay
                    section.data(100).logicalSrcIdx = 99;
                    section.data(100).dtTransOffset = 505;

                    ;% rtB.ndcqoyp0vx
                    section.data(101).logicalSrcIdx = 100;
                    section.data(101).dtTransOffset = 509;

                    ;% rtB.hv30h5tsxp
                    section.data(102).logicalSrcIdx = 101;
                    section.data(102).dtTransOffset = 512;

                    ;% rtB.iqlb144ibi
                    section.data(103).logicalSrcIdx = 102;
                    section.data(103).dtTransOffset = 513;

                    ;% rtB.pl32irduyc
                    section.data(104).logicalSrcIdx = 103;
                    section.data(104).dtTransOffset = 514;

                    ;% rtB.lgifpnih0e
                    section.data(105).logicalSrcIdx = 104;
                    section.data(105).dtTransOffset = 517;

                    ;% rtB.l0eiwonrmi
                    section.data(106).logicalSrcIdx = 105;
                    section.data(106).dtTransOffset = 520;

                    ;% rtB.hrtkouuc3i
                    section.data(107).logicalSrcIdx = 106;
                    section.data(107).dtTransOffset = 523;

                    ;% rtB.gmo4itlisv
                    section.data(108).logicalSrcIdx = 107;
                    section.data(108).dtTransOffset = 526;

                    ;% rtB.nexxpkogi1
                    section.data(109).logicalSrcIdx = 108;
                    section.data(109).dtTransOffset = 527;

                    ;% rtB.idwieepnv4
                    section.data(110).logicalSrcIdx = 109;
                    section.data(110).dtTransOffset = 530;

                    ;% rtB.h2vbgw30tf
                    section.data(111).logicalSrcIdx = 110;
                    section.data(111).dtTransOffset = 531;

                    ;% rtB.dndiwk2b1s
                    section.data(112).logicalSrcIdx = 111;
                    section.data(112).dtTransOffset = 534;

                    ;% rtB.gvavv3jbs0
                    section.data(113).logicalSrcIdx = 112;
                    section.data(113).dtTransOffset = 537;

                    ;% rtB.dw12uizmd2
                    section.data(114).logicalSrcIdx = 113;
                    section.data(114).dtTransOffset = 540;

                    ;% rtB.bdswxd3b20
                    section.data(115).logicalSrcIdx = 114;
                    section.data(115).dtTransOffset = 543;

                    ;% rtB.fh04srz5g0
                    section.data(116).logicalSrcIdx = 115;
                    section.data(116).dtTransOffset = 546;

                    ;% rtB.mkzi43umox
                    section.data(117).logicalSrcIdx = 116;
                    section.data(117).dtTransOffset = 547;

                    ;% rtB.cpzvyphms1
                    section.data(118).logicalSrcIdx = 117;
                    section.data(118).dtTransOffset = 550;

                    ;% rtB.oc1jrzdnhe
                    section.data(119).logicalSrcIdx = 118;
                    section.data(119).dtTransOffset = 559;

                    ;% rtB.hlnjhscgik
                    section.data(120).logicalSrcIdx = 119;
                    section.data(120).dtTransOffset = 560;

                    ;% rtB.fkzq14ohjb
                    section.data(121).logicalSrcIdx = 120;
                    section.data(121).dtTransOffset = 561;

                    ;% rtB.cafqksoz00
                    section.data(122).logicalSrcIdx = 121;
                    section.data(122).dtTransOffset = 565;

                    ;% rtB.c402o0md5n
                    section.data(123).logicalSrcIdx = 122;
                    section.data(123).dtTransOffset = 568;

                    ;% rtB.hwppgfjqwy
                    section.data(124).logicalSrcIdx = 123;
                    section.data(124).dtTransOffset = 571;

                    ;% rtB.i1lnvldtbv
                    section.data(125).logicalSrcIdx = 124;
                    section.data(125).dtTransOffset = 580;

                    ;% rtB.ahw1dydhgo
                    section.data(126).logicalSrcIdx = 125;
                    section.data(126).dtTransOffset = 589;

                    ;% rtB.ort4y5zlaf
                    section.data(127).logicalSrcIdx = 126;
                    section.data(127).dtTransOffset = 598;

                    ;% rtB.f1arg5nwa3
                    section.data(128).logicalSrcIdx = 127;
                    section.data(128).dtTransOffset = 601;

                    ;% rtB.eckq3h2tzb
                    section.data(129).logicalSrcIdx = 128;
                    section.data(129).dtTransOffset = 602;

                    ;% rtB.bsaqv0u04c
                    section.data(130).logicalSrcIdx = 129;
                    section.data(130).dtTransOffset = 605;

                    ;% rtB.alzuotumkg
                    section.data(131).logicalSrcIdx = 130;
                    section.data(131).dtTransOffset = 611;

                    ;% rtB.ley50ghvfj
                    section.data(132).logicalSrcIdx = 131;
                    section.data(132).dtTransOffset = 614;

                    ;% rtB.b0dwpzqyrm
                    section.data(133).logicalSrcIdx = 132;
                    section.data(133).dtTransOffset = 617;

                    ;% rtB.fbdf2tllgu
                    section.data(134).logicalSrcIdx = 133;
                    section.data(134).dtTransOffset = 620;

                    ;% rtB.my44vuapi0
                    section.data(135).logicalSrcIdx = 134;
                    section.data(135).dtTransOffset = 621;

                    ;% rtB.k3faf5atq0
                    section.data(136).logicalSrcIdx = 135;
                    section.data(136).dtTransOffset = 622;

                    ;% rtB.cg3eiq4i4s
                    section.data(137).logicalSrcIdx = 136;
                    section.data(137).dtTransOffset = 626;

                    ;% rtB.encqrxuu3y
                    section.data(138).logicalSrcIdx = 137;
                    section.data(138).dtTransOffset = 629;

                    ;% rtB.bdozk1fpny
                    section.data(139).logicalSrcIdx = 138;
                    section.data(139).dtTransOffset = 630;

                    ;% rtB.nbrrv2d3zg
                    section.data(140).logicalSrcIdx = 139;
                    section.data(140).dtTransOffset = 631;

                    ;% rtB.cuxk4j3qys
                    section.data(141).logicalSrcIdx = 140;
                    section.data(141).dtTransOffset = 634;

                    ;% rtB.eegu2ofo43
                    section.data(142).logicalSrcIdx = 141;
                    section.data(142).dtTransOffset = 637;

                    ;% rtB.hacw1qj3aj
                    section.data(143).logicalSrcIdx = 142;
                    section.data(143).dtTransOffset = 640;

                    ;% rtB.pxg5b2slaw
                    section.data(144).logicalSrcIdx = 143;
                    section.data(144).dtTransOffset = 643;

                    ;% rtB.h4h1ivsmcr
                    section.data(145).logicalSrcIdx = 144;
                    section.data(145).dtTransOffset = 644;

                    ;% rtB.jqzltcnumt
                    section.data(146).logicalSrcIdx = 145;
                    section.data(146).dtTransOffset = 647;

                    ;% rtB.etawd0whn4
                    section.data(147).logicalSrcIdx = 146;
                    section.data(147).dtTransOffset = 648;

                    ;% rtB.cbmxpbd2iu
                    section.data(148).logicalSrcIdx = 147;
                    section.data(148).dtTransOffset = 651;

                    ;% rtB.l1wqaovs4e
                    section.data(149).logicalSrcIdx = 148;
                    section.data(149).dtTransOffset = 654;

                    ;% rtB.f34tk41pk1
                    section.data(150).logicalSrcIdx = 149;
                    section.data(150).dtTransOffset = 657;

                    ;% rtB.kup2makpml
                    section.data(151).logicalSrcIdx = 150;
                    section.data(151).dtTransOffset = 660;

                    ;% rtB.d4it13dar5
                    section.data(152).logicalSrcIdx = 151;
                    section.data(152).dtTransOffset = 663;

                    ;% rtB.dl4tl2eodk
                    section.data(153).logicalSrcIdx = 152;
                    section.data(153).dtTransOffset = 664;

                    ;% rtB.obk5xjqd1g
                    section.data(154).logicalSrcIdx = 153;
                    section.data(154).dtTransOffset = 667;

                    ;% rtB.p00fanezjt
                    section.data(155).logicalSrcIdx = 154;
                    section.data(155).dtTransOffset = 676;

                    ;% rtB.b4b10et3s3
                    section.data(156).logicalSrcIdx = 155;
                    section.data(156).dtTransOffset = 677;

                    ;% rtB.flzta24bjr
                    section.data(157).logicalSrcIdx = 156;
                    section.data(157).dtTransOffset = 678;

                    ;% rtB.dv2svmenmy
                    section.data(158).logicalSrcIdx = 157;
                    section.data(158).dtTransOffset = 682;

                    ;% rtB.odlamrju54
                    section.data(159).logicalSrcIdx = 158;
                    section.data(159).dtTransOffset = 685;

                    ;% rtB.p3ru4cddgi
                    section.data(160).logicalSrcIdx = 159;
                    section.data(160).dtTransOffset = 688;

                    ;% rtB.bqjv3qidt0
                    section.data(161).logicalSrcIdx = 160;
                    section.data(161).dtTransOffset = 697;

                    ;% rtB.lkoo4fzueb
                    section.data(162).logicalSrcIdx = 161;
                    section.data(162).dtTransOffset = 706;

                    ;% rtB.gxha2ko2vz
                    section.data(163).logicalSrcIdx = 162;
                    section.data(163).dtTransOffset = 715;

                    ;% rtB.lffwcdbemm
                    section.data(164).logicalSrcIdx = 163;
                    section.data(164).dtTransOffset = 718;

                    ;% rtB.chig3jpwbn
                    section.data(165).logicalSrcIdx = 164;
                    section.data(165).dtTransOffset = 719;

                    ;% rtB.fcsl1vvaio
                    section.data(166).logicalSrcIdx = 165;
                    section.data(166).dtTransOffset = 722;

                    ;% rtB.pu15hoqemb
                    section.data(167).logicalSrcIdx = 166;
                    section.data(167).dtTransOffset = 728;

                    ;% rtB.emdpdwkzg3
                    section.data(168).logicalSrcIdx = 167;
                    section.data(168).dtTransOffset = 731;

                    ;% rtB.lxihdocxga
                    section.data(169).logicalSrcIdx = 168;
                    section.data(169).dtTransOffset = 734;

                    ;% rtB.nfak5a1vn4
                    section.data(170).logicalSrcIdx = 169;
                    section.data(170).dtTransOffset = 737;

                    ;% rtB.ntt11gc034
                    section.data(171).logicalSrcIdx = 170;
                    section.data(171).dtTransOffset = 738;

                    ;% rtB.hvd2sf2hag
                    section.data(172).logicalSrcIdx = 171;
                    section.data(172).dtTransOffset = 739;

                    ;% rtB.bpewjnssvt
                    section.data(173).logicalSrcIdx = 172;
                    section.data(173).dtTransOffset = 743;

                    ;% rtB.j1ydxe2usi
                    section.data(174).logicalSrcIdx = 173;
                    section.data(174).dtTransOffset = 746;

                    ;% rtB.ceb1rzrjkf
                    section.data(175).logicalSrcIdx = 174;
                    section.data(175).dtTransOffset = 747;

                    ;% rtB.lyhkcvfx5d
                    section.data(176).logicalSrcIdx = 175;
                    section.data(176).dtTransOffset = 748;

                    ;% rtB.enoycngmy3
                    section.data(177).logicalSrcIdx = 176;
                    section.data(177).dtTransOffset = 751;

                    ;% rtB.mw4hhpohry
                    section.data(178).logicalSrcIdx = 177;
                    section.data(178).dtTransOffset = 754;

                    ;% rtB.bhsuwz2d2l
                    section.data(179).logicalSrcIdx = 178;
                    section.data(179).dtTransOffset = 757;

                    ;% rtB.ezvzwoavxd
                    section.data(180).logicalSrcIdx = 179;
                    section.data(180).dtTransOffset = 760;

                    ;% rtB.msqcnwfh2s
                    section.data(181).logicalSrcIdx = 180;
                    section.data(181).dtTransOffset = 761;

                    ;% rtB.mp054zeo2w
                    section.data(182).logicalSrcIdx = 181;
                    section.data(182).dtTransOffset = 764;

                    ;% rtB.pfosbwya4i
                    section.data(183).logicalSrcIdx = 182;
                    section.data(183).dtTransOffset = 765;

                    ;% rtB.l01qgkxm4q
                    section.data(184).logicalSrcIdx = 183;
                    section.data(184).dtTransOffset = 768;

                    ;% rtB.babxk0cwue
                    section.data(185).logicalSrcIdx = 184;
                    section.data(185).dtTransOffset = 771;

                    ;% rtB.fzbgmnan4p
                    section.data(186).logicalSrcIdx = 185;
                    section.data(186).dtTransOffset = 774;

                    ;% rtB.oflmpiaez1
                    section.data(187).logicalSrcIdx = 186;
                    section.data(187).dtTransOffset = 777;

                    ;% rtB.efmfmwthqv
                    section.data(188).logicalSrcIdx = 187;
                    section.data(188).dtTransOffset = 780;

                    ;% rtB.ivaqfjshfj
                    section.data(189).logicalSrcIdx = 188;
                    section.data(189).dtTransOffset = 781;

                    ;% rtB.fzqka3pwrj
                    section.data(190).logicalSrcIdx = 189;
                    section.data(190).dtTransOffset = 784;

                    ;% rtB.egvr4ozujh
                    section.data(191).logicalSrcIdx = 190;
                    section.data(191).dtTransOffset = 793;

                    ;% rtB.cjrpkdpszw
                    section.data(192).logicalSrcIdx = 191;
                    section.data(192).dtTransOffset = 794;

                    ;% rtB.logpc2hca5
                    section.data(193).logicalSrcIdx = 192;
                    section.data(193).dtTransOffset = 795;

                    ;% rtB.ksrswfmi3x
                    section.data(194).logicalSrcIdx = 193;
                    section.data(194).dtTransOffset = 799;

                    ;% rtB.d2uj0bf0jq
                    section.data(195).logicalSrcIdx = 194;
                    section.data(195).dtTransOffset = 802;

                    ;% rtB.joimmmlrrj
                    section.data(196).logicalSrcIdx = 195;
                    section.data(196).dtTransOffset = 805;

                    ;% rtB.aaw4gtn4yc
                    section.data(197).logicalSrcIdx = 196;
                    section.data(197).dtTransOffset = 814;

                    ;% rtB.jnz14mmcr3
                    section.data(198).logicalSrcIdx = 197;
                    section.data(198).dtTransOffset = 823;

                    ;% rtB.le0lvwzojd
                    section.data(199).logicalSrcIdx = 198;
                    section.data(199).dtTransOffset = 832;

                    ;% rtB.k3vxq1m344
                    section.data(200).logicalSrcIdx = 199;
                    section.data(200).dtTransOffset = 835;

                    ;% rtB.e3hmodlujb
                    section.data(201).logicalSrcIdx = 200;
                    section.data(201).dtTransOffset = 836;

                    ;% rtB.b3yujzbge2
                    section.data(202).logicalSrcIdx = 201;
                    section.data(202).dtTransOffset = 839;

                    ;% rtB.gf5n1xgexk
                    section.data(203).logicalSrcIdx = 202;
                    section.data(203).dtTransOffset = 845;

                    ;% rtB.e3buvmif1c
                    section.data(204).logicalSrcIdx = 203;
                    section.data(204).dtTransOffset = 848;

                    ;% rtB.j4fldanieq
                    section.data(205).logicalSrcIdx = 204;
                    section.data(205).dtTransOffset = 851;

                    ;% rtB.cdgxtvlobx
                    section.data(206).logicalSrcIdx = 205;
                    section.data(206).dtTransOffset = 854;

                    ;% rtB.lbmtuwi01h
                    section.data(207).logicalSrcIdx = 206;
                    section.data(207).dtTransOffset = 855;

                    ;% rtB.ieyrepmb05
                    section.data(208).logicalSrcIdx = 207;
                    section.data(208).dtTransOffset = 856;

                    ;% rtB.lyjnnzstpa
                    section.data(209).logicalSrcIdx = 208;
                    section.data(209).dtTransOffset = 857;

                    ;% rtB.favu2vq25u
                    section.data(210).logicalSrcIdx = 209;
                    section.data(210).dtTransOffset = 861;

                    ;% rtB.h1zdr1aani
                    section.data(211).logicalSrcIdx = 210;
                    section.data(211).dtTransOffset = 864;

                    ;% rtB.m4bi5zptkh
                    section.data(212).logicalSrcIdx = 211;
                    section.data(212).dtTransOffset = 867;

                    ;% rtB.cab5adcq1r
                    section.data(213).logicalSrcIdx = 212;
                    section.data(213).dtTransOffset = 870;

                    ;% rtB.avcgwoxwhd
                    section.data(214).logicalSrcIdx = 213;
                    section.data(214).dtTransOffset = 871;

                    ;% rtB.f4nnctaymu
                    section.data(215).logicalSrcIdx = 214;
                    section.data(215).dtTransOffset = 872;

                    ;% rtB.mcasp4wuog
                    section.data(216).logicalSrcIdx = 215;
                    section.data(216).dtTransOffset = 875;

                    ;% rtB.mrdivqgmvm
                    section.data(217).logicalSrcIdx = 216;
                    section.data(217).dtTransOffset = 878;

                    ;% rtB.g3vpjafupq
                    section.data(218).logicalSrcIdx = 217;
                    section.data(218).dtTransOffset = 881;

                    ;% rtB.btylsma43j
                    section.data(219).logicalSrcIdx = 218;
                    section.data(219).dtTransOffset = 884;

                    ;% rtB.my424itvki
                    section.data(220).logicalSrcIdx = 219;
                    section.data(220).dtTransOffset = 885;

                    ;% rtB.fmqgelgm0i
                    section.data(221).logicalSrcIdx = 220;
                    section.data(221).dtTransOffset = 888;

                    ;% rtB.gmivsxqevf
                    section.data(222).logicalSrcIdx = 221;
                    section.data(222).dtTransOffset = 891;

                    ;% rtB.cxwqk2amda
                    section.data(223).logicalSrcIdx = 222;
                    section.data(223).dtTransOffset = 894;

                    ;% rtB.jgkbpstejl
                    section.data(224).logicalSrcIdx = 223;
                    section.data(224).dtTransOffset = 897;

                    ;% rtB.evftq3sk0q
                    section.data(225).logicalSrcIdx = 224;
                    section.data(225).dtTransOffset = 898;

                    ;% rtB.hnzvyzmoy4
                    section.data(226).logicalSrcIdx = 225;
                    section.data(226).dtTransOffset = 901;

                    ;% rtB.gtu1mufolx
                    section.data(227).logicalSrcIdx = 226;
                    section.data(227).dtTransOffset = 910;

                    ;% rtB.ifl2ongx1d
                    section.data(228).logicalSrcIdx = 227;
                    section.data(228).dtTransOffset = 911;

                    ;% rtB.a3qlbfbxxo
                    section.data(229).logicalSrcIdx = 228;
                    section.data(229).dtTransOffset = 912;

                    ;% rtB.glnk1qsnt1
                    section.data(230).logicalSrcIdx = 229;
                    section.data(230).dtTransOffset = 916;

                    ;% rtB.li0pbutf3d
                    section.data(231).logicalSrcIdx = 230;
                    section.data(231).dtTransOffset = 919;

                    ;% rtB.nzexd0tfyn
                    section.data(232).logicalSrcIdx = 231;
                    section.data(232).dtTransOffset = 922;

                    ;% rtB.lrkvtt341w
                    section.data(233).logicalSrcIdx = 232;
                    section.data(233).dtTransOffset = 931;

                    ;% rtB.kptzapguuy
                    section.data(234).logicalSrcIdx = 233;
                    section.data(234).dtTransOffset = 940;

                    ;% rtB.nwco0rdsmq
                    section.data(235).logicalSrcIdx = 234;
                    section.data(235).dtTransOffset = 949;

                    ;% rtB.j2h4dmgm1o
                    section.data(236).logicalSrcIdx = 235;
                    section.data(236).dtTransOffset = 952;

                    ;% rtB.kdzk5onxbh
                    section.data(237).logicalSrcIdx = 236;
                    section.data(237).dtTransOffset = 953;

                    ;% rtB.oiy31a51av
                    section.data(238).logicalSrcIdx = 237;
                    section.data(238).dtTransOffset = 956;

                    ;% rtB.athb42mz2d
                    section.data(239).logicalSrcIdx = 238;
                    section.data(239).dtTransOffset = 962;

                    ;% rtB.exqp2uohe5
                    section.data(240).logicalSrcIdx = 239;
                    section.data(240).dtTransOffset = 965;

                    ;% rtB.o2znllyvoq
                    section.data(241).logicalSrcIdx = 240;
                    section.data(241).dtTransOffset = 968;

                    ;% rtB.bbtvepo4p4
                    section.data(242).logicalSrcIdx = 241;
                    section.data(242).dtTransOffset = 971;

                    ;% rtB.kq5inrtfnq
                    section.data(243).logicalSrcIdx = 242;
                    section.data(243).dtTransOffset = 972;

                    ;% rtB.c00wpmk13d
                    section.data(244).logicalSrcIdx = 243;
                    section.data(244).dtTransOffset = 973;

                    ;% rtB.huacygyj0o
                    section.data(245).logicalSrcIdx = 244;
                    section.data(245).dtTransOffset = 974;

                    ;% rtB.gg2oxfw35t
                    section.data(246).logicalSrcIdx = 245;
                    section.data(246).dtTransOffset = 978;

                    ;% rtB.lperucgciu
                    section.data(247).logicalSrcIdx = 246;
                    section.data(247).dtTransOffset = 981;

                    ;% rtB.jqxy4sgz1r
                    section.data(248).logicalSrcIdx = 247;
                    section.data(248).dtTransOffset = 984;

                    ;% rtB.aizglk3yjb
                    section.data(249).logicalSrcIdx = 248;
                    section.data(249).dtTransOffset = 987;

                    ;% rtB.avjauevoe4
                    section.data(250).logicalSrcIdx = 249;
                    section.data(250).dtTransOffset = 988;

                    ;% rtB.azxdea4e00
                    section.data(251).logicalSrcIdx = 250;
                    section.data(251).dtTransOffset = 989;

                    ;% rtB.afxh4mt3gu
                    section.data(252).logicalSrcIdx = 251;
                    section.data(252).dtTransOffset = 992;

                    ;% rtB.gdhtjuaiqp
                    section.data(253).logicalSrcIdx = 252;
                    section.data(253).dtTransOffset = 995;

                    ;% rtB.e0gxuicxzo
                    section.data(254).logicalSrcIdx = 253;
                    section.data(254).dtTransOffset = 998;

                    ;% rtB.kbih5351xy
                    section.data(255).logicalSrcIdx = 254;
                    section.data(255).dtTransOffset = 999;

                    ;% rtB.j2tnneegmr
                    section.data(256).logicalSrcIdx = 255;
                    section.data(256).dtTransOffset = 1002;

                    ;% rtB.c513o5pg5r
                    section.data(257).logicalSrcIdx = 256;
                    section.data(257).dtTransOffset = 1005;

                    ;% rtB.pfmi5qvzhf
                    section.data(258).logicalSrcIdx = 257;
                    section.data(258).dtTransOffset = 1009;

                    ;% rtB.gi35bo1nw4
                    section.data(259).logicalSrcIdx = 258;
                    section.data(259).dtTransOffset = 1013;

                    ;% rtB.iii5xyl5wn
                    section.data(260).logicalSrcIdx = 259;
                    section.data(260).dtTransOffset = 1017;

                    ;% rtB.mjqwg1vruy
                    section.data(261).logicalSrcIdx = 260;
                    section.data(261).dtTransOffset = 1021;

                    ;% rtB.eo0yhyzcnn
                    section.data(262).logicalSrcIdx = 261;
                    section.data(262).dtTransOffset = 1025;

                    ;% rtB.c5dvbg4axm
                    section.data(263).logicalSrcIdx = 262;
                    section.data(263).dtTransOffset = 1029;

                    ;% rtB.bxbu5xamaf
                    section.data(264).logicalSrcIdx = 263;
                    section.data(264).dtTransOffset = 1033;

                    ;% rtB.foqtyz1zpt
                    section.data(265).logicalSrcIdx = 264;
                    section.data(265).dtTransOffset = 1037;

                    ;% rtB.gy3u0julyl
                    section.data(266).logicalSrcIdx = 265;
                    section.data(266).dtTransOffset = 1041;

                    ;% rtB.hytxiurw0s
                    section.data(267).logicalSrcIdx = 266;
                    section.data(267).dtTransOffset = 1045;

                    ;% rtB.nubcraskrw
                    section.data(268).logicalSrcIdx = 267;
                    section.data(268).dtTransOffset = 1049;

                    ;% rtB.lp0cxari03
                    section.data(269).logicalSrcIdx = 268;
                    section.data(269).dtTransOffset = 1053;

                    ;% rtB.phx5lsta5w
                    section.data(270).logicalSrcIdx = 269;
                    section.data(270).dtTransOffset = 1057;

                    ;% rtB.hoap0zr3eu
                    section.data(271).logicalSrcIdx = 270;
                    section.data(271).dtTransOffset = 1061;

                    ;% rtB.enp5f3s002
                    section.data(272).logicalSrcIdx = 271;
                    section.data(272).dtTransOffset = 1065;

                    ;% rtB.nsl31ekm0f
                    section.data(273).logicalSrcIdx = 272;
                    section.data(273).dtTransOffset = 1069;

                    ;% rtB.hrqvwki2jl
                    section.data(274).logicalSrcIdx = 273;
                    section.data(274).dtTransOffset = 1073;

                    ;% rtB.dypcby02yn
                    section.data(275).logicalSrcIdx = 274;
                    section.data(275).dtTransOffset = 1077;

                    ;% rtB.a100glwqkz
                    section.data(276).logicalSrcIdx = 275;
                    section.data(276).dtTransOffset = 1081;

                    ;% rtB.c3twv1nzdj
                    section.data(277).logicalSrcIdx = 276;
                    section.data(277).dtTransOffset = 1085;

                    ;% rtB.iurh4h1cla
                    section.data(278).logicalSrcIdx = 277;
                    section.data(278).dtTransOffset = 1089;

                    ;% rtB.har0a2wrmi
                    section.data(279).logicalSrcIdx = 278;
                    section.data(279).dtTransOffset = 1093;

                    ;% rtB.ib2ieuf2an
                    section.data(280).logicalSrcIdx = 279;
                    section.data(280).dtTransOffset = 1097;

                    ;% rtB.ojpovrodvr
                    section.data(281).logicalSrcIdx = 280;
                    section.data(281).dtTransOffset = 1101;

                    ;% rtB.fcgmpow5uz
                    section.data(282).logicalSrcIdx = 281;
                    section.data(282).dtTransOffset = 1105;

                    ;% rtB.jylntbxfpj
                    section.data(283).logicalSrcIdx = 282;
                    section.data(283).dtTransOffset = 1109;

                    ;% rtB.g1yo2p3pxt
                    section.data(284).logicalSrcIdx = 283;
                    section.data(284).dtTransOffset = 1113;

                    ;% rtB.az4eo3zh5o
                    section.data(285).logicalSrcIdx = 284;
                    section.data(285).dtTransOffset = 1117;

                    ;% rtB.ch1tbnhyon
                    section.data(286).logicalSrcIdx = 285;
                    section.data(286).dtTransOffset = 1121;

                    ;% rtB.eoseoj2f3b
                    section.data(287).logicalSrcIdx = 286;
                    section.data(287).dtTransOffset = 1125;

                    ;% rtB.aigpe1li2p
                    section.data(288).logicalSrcIdx = 287;
                    section.data(288).dtTransOffset = 1129;

                    ;% rtB.pbn2nzk2k3
                    section.data(289).logicalSrcIdx = 288;
                    section.data(289).dtTransOffset = 1133;

                    ;% rtB.hffhwq1imf
                    section.data(290).logicalSrcIdx = 289;
                    section.data(290).dtTransOffset = 1137;

                    ;% rtB.gbxuu4wtte
                    section.data(291).logicalSrcIdx = 290;
                    section.data(291).dtTransOffset = 1141;

                    ;% rtB.k5nlme3kns
                    section.data(292).logicalSrcIdx = 291;
                    section.data(292).dtTransOffset = 1145;

                    ;% rtB.j52cfj2d1w
                    section.data(293).logicalSrcIdx = 292;
                    section.data(293).dtTransOffset = 1149;

                    ;% rtB.ftvz2ekipg
                    section.data(294).logicalSrcIdx = 293;
                    section.data(294).dtTransOffset = 1153;

                    ;% rtB.awgskjjmsb
                    section.data(295).logicalSrcIdx = 294;
                    section.data(295).dtTransOffset = 1157;

                    ;% rtB.fl4rvtkarn
                    section.data(296).logicalSrcIdx = 295;
                    section.data(296).dtTransOffset = 1161;

                    ;% rtB.jp22vgjwwm
                    section.data(297).logicalSrcIdx = 296;
                    section.data(297).dtTransOffset = 1165;

                    ;% rtB.pzpjv5s5n2
                    section.data(298).logicalSrcIdx = 297;
                    section.data(298).dtTransOffset = 1169;

                    ;% rtB.nsmqj1yv3h
                    section.data(299).logicalSrcIdx = 298;
                    section.data(299).dtTransOffset = 1173;

                    ;% rtB.mdrxplc5qx
                    section.data(300).logicalSrcIdx = 299;
                    section.data(300).dtTransOffset = 1176;

                    ;% rtB.idcqc0uibn
                    section.data(301).logicalSrcIdx = 300;
                    section.data(301).dtTransOffset = 1179;

                    ;% rtB.lkwtwvkutl
                    section.data(302).logicalSrcIdx = 301;
                    section.data(302).dtTransOffset = 1182;

                    ;% rtB.dunjuwexfg
                    section.data(303).logicalSrcIdx = 302;
                    section.data(303).dtTransOffset = 1185;

                    ;% rtB.oimhqzeh2t
                    section.data(304).logicalSrcIdx = 303;
                    section.data(304).dtTransOffset = 1188;

                    ;% rtB.bcudxny0c0
                    section.data(305).logicalSrcIdx = 304;
                    section.data(305).dtTransOffset = 1191;

                    ;% rtB.ibott2adsg
                    section.data(306).logicalSrcIdx = 305;
                    section.data(306).dtTransOffset = 1194;

                    ;% rtB.ovchpiv2ws
                    section.data(307).logicalSrcIdx = 306;
                    section.data(307).dtTransOffset = 1197;

                    ;% rtB.eejq0hbv30
                    section.data(308).logicalSrcIdx = 307;
                    section.data(308).dtTransOffset = 1200;

                    ;% rtB.joumqp0zpt
                    section.data(309).logicalSrcIdx = 308;
                    section.data(309).dtTransOffset = 1203;

                    ;% rtB.l02fn5ndrj
                    section.data(310).logicalSrcIdx = 309;
                    section.data(310).dtTransOffset = 1206;

                    ;% rtB.emllridil0
                    section.data(311).logicalSrcIdx = 310;
                    section.data(311).dtTransOffset = 1209;

                    ;% rtB.d1zqfcz4zo
                    section.data(312).logicalSrcIdx = 311;
                    section.data(312).dtTransOffset = 1212;

                    ;% rtB.bjhxaihh05
                    section.data(313).logicalSrcIdx = 312;
                    section.data(313).dtTransOffset = 1215;

                    ;% rtB.jltrdwyzmw
                    section.data(314).logicalSrcIdx = 313;
                    section.data(314).dtTransOffset = 1218;

                    ;% rtB.fc15cwanga
                    section.data(315).logicalSrcIdx = 314;
                    section.data(315).dtTransOffset = 1221;

                    ;% rtB.lebemiabnp
                    section.data(316).logicalSrcIdx = 315;
                    section.data(316).dtTransOffset = 1224;

                    ;% rtB.p5ihfgrspp
                    section.data(317).logicalSrcIdx = 316;
                    section.data(317).dtTransOffset = 1227;

                    ;% rtB.itx5qkpva4
                    section.data(318).logicalSrcIdx = 317;
                    section.data(318).dtTransOffset = 1230;

                    ;% rtB.e4tk4xaa2e
                    section.data(319).logicalSrcIdx = 318;
                    section.data(319).dtTransOffset = 1233;

                    ;% rtB.mdxnlbsfvm
                    section.data(320).logicalSrcIdx = 319;
                    section.data(320).dtTransOffset = 1236;

                    ;% rtB.fvtfdqea2r
                    section.data(321).logicalSrcIdx = 320;
                    section.data(321).dtTransOffset = 1239;

                    ;% rtB.inhs5zftab
                    section.data(322).logicalSrcIdx = 321;
                    section.data(322).dtTransOffset = 1242;

                    ;% rtB.hrudfjbz3i
                    section.data(323).logicalSrcIdx = 322;
                    section.data(323).dtTransOffset = 1245;

                    ;% rtB.cschdiniut
                    section.data(324).logicalSrcIdx = 323;
                    section.data(324).dtTransOffset = 1248;

                    ;% rtB.o4wni3515w
                    section.data(325).logicalSrcIdx = 324;
                    section.data(325).dtTransOffset = 1251;

                    ;% rtB.p5hpxiz5u3
                    section.data(326).logicalSrcIdx = 325;
                    section.data(326).dtTransOffset = 1254;

                    ;% rtB.h2ica3jpzg
                    section.data(327).logicalSrcIdx = 326;
                    section.data(327).dtTransOffset = 1257;

                    ;% rtB.jqlckttaam
                    section.data(328).logicalSrcIdx = 327;
                    section.data(328).dtTransOffset = 1260;

                    ;% rtB.npob1uuqas
                    section.data(329).logicalSrcIdx = 328;
                    section.data(329).dtTransOffset = 1263;

                    ;% rtB.itxt4ytwah
                    section.data(330).logicalSrcIdx = 329;
                    section.data(330).dtTransOffset = 1266;

                    ;% rtB.nl3crkyc5i
                    section.data(331).logicalSrcIdx = 330;
                    section.data(331).dtTransOffset = 1269;

                    ;% rtB.kdrism1ogn
                    section.data(332).logicalSrcIdx = 331;
                    section.data(332).dtTransOffset = 1272;

                    ;% rtB.nuca2ykly0
                    section.data(333).logicalSrcIdx = 332;
                    section.data(333).dtTransOffset = 1275;

            nTotData = nTotData + section.nData;
            sigMap.sections(1) = section;
            clear section


            ;%
            ;% Non-auto Data (signal)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        sigMap.nTotData = nTotData;



    ;%*******************
    ;% Create DWork Map *
    ;%*******************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 5;
        sectIdxOffset = 1;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc dworkMap
        ;%
        dworkMap.nSections           = nTotSects;
        dworkMap.sectIdxOffset       = sectIdxOffset;
            dworkMap.sections(nTotSects) = dumSection; %prealloc
        dworkMap.nTotData            = -1;

        ;%
        ;% Auto data (rtDW)
        ;%
            section.nData     = 44;
            section.data(44)  = dumData; %prealloc

                    ;% rtDW.ik1trwr3gq
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.m0rr22lxax
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 2;

                    ;% rtDW.jhm0cqabg3
                    section.data(3).logicalSrcIdx = 2;
                    section.data(3).dtTransOffset = 4;

                    ;% rtDW.oq35dv23vw
                    section.data(4).logicalSrcIdx = 3;
                    section.data(4).dtTransOffset = 6;

                    ;% rtDW.lwhmmfv2ps
                    section.data(5).logicalSrcIdx = 4;
                    section.data(5).dtTransOffset = 8;

                    ;% rtDW.iyivz4czyo
                    section.data(6).logicalSrcIdx = 5;
                    section.data(6).dtTransOffset = 10;

                    ;% rtDW.ls5rwnwukt
                    section.data(7).logicalSrcIdx = 6;
                    section.data(7).dtTransOffset = 12;

                    ;% rtDW.mebb30zlfb
                    section.data(8).logicalSrcIdx = 7;
                    section.data(8).dtTransOffset = 14;

                    ;% rtDW.kmluc3bnxm
                    section.data(9).logicalSrcIdx = 8;
                    section.data(9).dtTransOffset = 16;

                    ;% rtDW.ardpmlcavq
                    section.data(10).logicalSrcIdx = 9;
                    section.data(10).dtTransOffset = 18;

                    ;% rtDW.o2te45t5iz
                    section.data(11).logicalSrcIdx = 10;
                    section.data(11).dtTransOffset = 20;

                    ;% rtDW.ceunniigow
                    section.data(12).logicalSrcIdx = 11;
                    section.data(12).dtTransOffset = 22;

                    ;% rtDW.by4chuqxsd
                    section.data(13).logicalSrcIdx = 12;
                    section.data(13).dtTransOffset = 24;

                    ;% rtDW.pjygrikwyo
                    section.data(14).logicalSrcIdx = 13;
                    section.data(14).dtTransOffset = 26;

                    ;% rtDW.dbasyf4pfi
                    section.data(15).logicalSrcIdx = 14;
                    section.data(15).dtTransOffset = 28;

                    ;% rtDW.nj2tur0kwz
                    section.data(16).logicalSrcIdx = 15;
                    section.data(16).dtTransOffset = 30;

                    ;% rtDW.f1cbvtmv22
                    section.data(17).logicalSrcIdx = 16;
                    section.data(17).dtTransOffset = 32;

                    ;% rtDW.mh3gjthtbl
                    section.data(18).logicalSrcIdx = 17;
                    section.data(18).dtTransOffset = 34;

                    ;% rtDW.jwgslfbtic
                    section.data(19).logicalSrcIdx = 18;
                    section.data(19).dtTransOffset = 36;

                    ;% rtDW.ezstggpf2j
                    section.data(20).logicalSrcIdx = 19;
                    section.data(20).dtTransOffset = 38;

                    ;% rtDW.pd45haa5mh
                    section.data(21).logicalSrcIdx = 20;
                    section.data(21).dtTransOffset = 40;

                    ;% rtDW.lngqcivqqg
                    section.data(22).logicalSrcIdx = 21;
                    section.data(22).dtTransOffset = 42;

                    ;% rtDW.ozmcvi5b31
                    section.data(23).logicalSrcIdx = 22;
                    section.data(23).dtTransOffset = 44;

                    ;% rtDW.hdadv1kifl
                    section.data(24).logicalSrcIdx = 23;
                    section.data(24).dtTransOffset = 46;

                    ;% rtDW.dhds3swxi0
                    section.data(25).logicalSrcIdx = 24;
                    section.data(25).dtTransOffset = 48;

                    ;% rtDW.oqc0b5does
                    section.data(26).logicalSrcIdx = 25;
                    section.data(26).dtTransOffset = 50;

                    ;% rtDW.ieohtob00g
                    section.data(27).logicalSrcIdx = 26;
                    section.data(27).dtTransOffset = 52;

                    ;% rtDW.eiyeel4sy0
                    section.data(28).logicalSrcIdx = 27;
                    section.data(28).dtTransOffset = 54;

                    ;% rtDW.lgu0i43xpl
                    section.data(29).logicalSrcIdx = 28;
                    section.data(29).dtTransOffset = 56;

                    ;% rtDW.n01cjpi1wr
                    section.data(30).logicalSrcIdx = 29;
                    section.data(30).dtTransOffset = 58;

                    ;% rtDW.bviufihvoa
                    section.data(31).logicalSrcIdx = 30;
                    section.data(31).dtTransOffset = 60;

                    ;% rtDW.i0bfhtiky2
                    section.data(32).logicalSrcIdx = 31;
                    section.data(32).dtTransOffset = 62;

                    ;% rtDW.lq23ppo5t1
                    section.data(33).logicalSrcIdx = 32;
                    section.data(33).dtTransOffset = 64;

                    ;% rtDW.ffetch5fn3
                    section.data(34).logicalSrcIdx = 33;
                    section.data(34).dtTransOffset = 66;

                    ;% rtDW.aj3ho3s0dk
                    section.data(35).logicalSrcIdx = 34;
                    section.data(35).dtTransOffset = 68;

                    ;% rtDW.jjxvwteapb
                    section.data(36).logicalSrcIdx = 35;
                    section.data(36).dtTransOffset = 70;

                    ;% rtDW.mqctzxdg0o
                    section.data(37).logicalSrcIdx = 36;
                    section.data(37).dtTransOffset = 72;

                    ;% rtDW.nnmjponkr5
                    section.data(38).logicalSrcIdx = 37;
                    section.data(38).dtTransOffset = 74;

                    ;% rtDW.ky3lis0a5h
                    section.data(39).logicalSrcIdx = 38;
                    section.data(39).dtTransOffset = 76;

                    ;% rtDW.iftqygn0yl
                    section.data(40).logicalSrcIdx = 39;
                    section.data(40).dtTransOffset = 78;

                    ;% rtDW.ndljvxaegf
                    section.data(41).logicalSrcIdx = 40;
                    section.data(41).dtTransOffset = 80;

                    ;% rtDW.k5hl2nqzw5
                    section.data(42).logicalSrcIdx = 41;
                    section.data(42).dtTransOffset = 82;

                    ;% rtDW.exk4rrasy3
                    section.data(43).logicalSrcIdx = 42;
                    section.data(43).dtTransOffset = 84;

                    ;% rtDW.oz0f0tgc3u
                    section.data(44).logicalSrcIdx = 43;
                    section.data(44).dtTransOffset = 85;

            nTotData = nTotData + section.nData;
            dworkMap.sections(1) = section;
            clear section

            section.nData     = 14;
            section.data(14)  = dumData; %prealloc

                    ;% rtDW.ctycytv1sy
                    section.data(1).logicalSrcIdx = 44;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.m2wgpmvjvy
                    section.data(2).logicalSrcIdx = 45;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.oqiguspxfn
                    section.data(3).logicalSrcIdx = 46;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.bigiskgqfj
                    section.data(4).logicalSrcIdx = 47;
                    section.data(4).dtTransOffset = 3;

                    ;% rtDW.djxei35le3
                    section.data(5).logicalSrcIdx = 48;
                    section.data(5).dtTransOffset = 4;

                    ;% rtDW.dikp21y4fz
                    section.data(6).logicalSrcIdx = 49;
                    section.data(6).dtTransOffset = 5;

                    ;% rtDW.fqoysoyi51
                    section.data(7).logicalSrcIdx = 50;
                    section.data(7).dtTransOffset = 6;

                    ;% rtDW.dcyq5ws3qs
                    section.data(8).logicalSrcIdx = 51;
                    section.data(8).dtTransOffset = 7;

                    ;% rtDW.ixn514xqve
                    section.data(9).logicalSrcIdx = 52;
                    section.data(9).dtTransOffset = 8;

                    ;% rtDW.baw0vps5za
                    section.data(10).logicalSrcIdx = 53;
                    section.data(10).dtTransOffset = 9;

                    ;% rtDW.l4r0r31ubh.AQHandles
                    section.data(11).logicalSrcIdx = 54;
                    section.data(11).dtTransOffset = 10;

                    ;% rtDW.bwezjkpzim
                    section.data(12).logicalSrcIdx = 55;
                    section.data(12).dtTransOffset = 11;

                    ;% rtDW.bi5anfwza5
                    section.data(13).logicalSrcIdx = 56;
                    section.data(13).dtTransOffset = 12;

                    ;% rtDW.gyisyfvd4m
                    section.data(14).logicalSrcIdx = 57;
                    section.data(14).dtTransOffset = 13;

            nTotData = nTotData + section.nData;
            dworkMap.sections(2) = section;
            clear section

            section.nData     = 2;
            section.data(2)  = dumData; %prealloc

                    ;% rtDW.h2bg1hotvh
                    section.data(1).logicalSrcIdx = 58;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.glrxfbh4lh
                    section.data(2).logicalSrcIdx = 59;
                    section.data(2).dtTransOffset = 1;

            nTotData = nTotData + section.nData;
            dworkMap.sections(3) = section;
            clear section

            section.nData     = 42;
            section.data(42)  = dumData; %prealloc

                    ;% rtDW.oay4r2t35k
                    section.data(1).logicalSrcIdx = 60;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.itbalhm53p
                    section.data(2).logicalSrcIdx = 61;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.o3hz2jqijw
                    section.data(3).logicalSrcIdx = 62;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.a1ctm05tly
                    section.data(4).logicalSrcIdx = 63;
                    section.data(4).dtTransOffset = 3;

                    ;% rtDW.gib2w1lke4
                    section.data(5).logicalSrcIdx = 64;
                    section.data(5).dtTransOffset = 4;

                    ;% rtDW.ncvj32mtr4
                    section.data(6).logicalSrcIdx = 65;
                    section.data(6).dtTransOffset = 5;

                    ;% rtDW.ojwdvmmbzv
                    section.data(7).logicalSrcIdx = 66;
                    section.data(7).dtTransOffset = 6;

                    ;% rtDW.blmbxswg4c
                    section.data(8).logicalSrcIdx = 67;
                    section.data(8).dtTransOffset = 7;

                    ;% rtDW.eiqm0ndfyg
                    section.data(9).logicalSrcIdx = 68;
                    section.data(9).dtTransOffset = 8;

                    ;% rtDW.pibyff0eql
                    section.data(10).logicalSrcIdx = 69;
                    section.data(10).dtTransOffset = 9;

                    ;% rtDW.l2ec5jammp
                    section.data(11).logicalSrcIdx = 70;
                    section.data(11).dtTransOffset = 10;

                    ;% rtDW.iuadxosr40
                    section.data(12).logicalSrcIdx = 71;
                    section.data(12).dtTransOffset = 11;

                    ;% rtDW.iqaxcni1cl
                    section.data(13).logicalSrcIdx = 72;
                    section.data(13).dtTransOffset = 12;

                    ;% rtDW.o4w412lyhz
                    section.data(14).logicalSrcIdx = 73;
                    section.data(14).dtTransOffset = 13;

                    ;% rtDW.l4r1nqhb1g
                    section.data(15).logicalSrcIdx = 74;
                    section.data(15).dtTransOffset = 14;

                    ;% rtDW.gq3uhkidnv
                    section.data(16).logicalSrcIdx = 75;
                    section.data(16).dtTransOffset = 15;

                    ;% rtDW.emuswyjmhs
                    section.data(17).logicalSrcIdx = 76;
                    section.data(17).dtTransOffset = 16;

                    ;% rtDW.eqmbwqv4uu
                    section.data(18).logicalSrcIdx = 77;
                    section.data(18).dtTransOffset = 17;

                    ;% rtDW.bxfixvxfdw
                    section.data(19).logicalSrcIdx = 78;
                    section.data(19).dtTransOffset = 18;

                    ;% rtDW.dei1252yrl
                    section.data(20).logicalSrcIdx = 79;
                    section.data(20).dtTransOffset = 19;

                    ;% rtDW.ogqlv3h0ke
                    section.data(21).logicalSrcIdx = 80;
                    section.data(21).dtTransOffset = 20;

                    ;% rtDW.cpgojctm31
                    section.data(22).logicalSrcIdx = 81;
                    section.data(22).dtTransOffset = 21;

                    ;% rtDW.hgiqtwygy1
                    section.data(23).logicalSrcIdx = 82;
                    section.data(23).dtTransOffset = 22;

                    ;% rtDW.ddu0m0coos
                    section.data(24).logicalSrcIdx = 83;
                    section.data(24).dtTransOffset = 23;

                    ;% rtDW.plosltfnse
                    section.data(25).logicalSrcIdx = 84;
                    section.data(25).dtTransOffset = 24;

                    ;% rtDW.i2wzmz255s
                    section.data(26).logicalSrcIdx = 85;
                    section.data(26).dtTransOffset = 25;

                    ;% rtDW.a1er4qtcuo
                    section.data(27).logicalSrcIdx = 86;
                    section.data(27).dtTransOffset = 26;

                    ;% rtDW.chs3ptplyf
                    section.data(28).logicalSrcIdx = 87;
                    section.data(28).dtTransOffset = 27;

                    ;% rtDW.jtymq0r3wm
                    section.data(29).logicalSrcIdx = 88;
                    section.data(29).dtTransOffset = 28;

                    ;% rtDW.avqe5nr4ga
                    section.data(30).logicalSrcIdx = 89;
                    section.data(30).dtTransOffset = 29;

                    ;% rtDW.h53r4zk0cj
                    section.data(31).logicalSrcIdx = 90;
                    section.data(31).dtTransOffset = 30;

                    ;% rtDW.gktfawr2st
                    section.data(32).logicalSrcIdx = 91;
                    section.data(32).dtTransOffset = 31;

                    ;% rtDW.bj4ejpvtiv
                    section.data(33).logicalSrcIdx = 92;
                    section.data(33).dtTransOffset = 32;

                    ;% rtDW.eyhr14ycpb
                    section.data(34).logicalSrcIdx = 93;
                    section.data(34).dtTransOffset = 33;

                    ;% rtDW.i0q150xsvv
                    section.data(35).logicalSrcIdx = 94;
                    section.data(35).dtTransOffset = 34;

                    ;% rtDW.inbjvmxdis
                    section.data(36).logicalSrcIdx = 95;
                    section.data(36).dtTransOffset = 35;

                    ;% rtDW.lam1o4zpo4
                    section.data(37).logicalSrcIdx = 96;
                    section.data(37).dtTransOffset = 36;

                    ;% rtDW.ot4rduhsnb
                    section.data(38).logicalSrcIdx = 97;
                    section.data(38).dtTransOffset = 37;

                    ;% rtDW.gmlynrr0kl
                    section.data(39).logicalSrcIdx = 98;
                    section.data(39).dtTransOffset = 38;

                    ;% rtDW.dfembtpulp
                    section.data(40).logicalSrcIdx = 99;
                    section.data(40).dtTransOffset = 39;

                    ;% rtDW.mg55yfehmm
                    section.data(41).logicalSrcIdx = 100;
                    section.data(41).dtTransOffset = 40;

                    ;% rtDW.knhfhd3fs2
                    section.data(42).logicalSrcIdx = 101;
                    section.data(42).dtTransOffset = 41;

            nTotData = nTotData + section.nData;
            dworkMap.sections(4) = section;
            clear section

            section.nData     = 30;
            section.data(30)  = dumData; %prealloc

                    ;% rtDW.l4gngyaono
                    section.data(1).logicalSrcIdx = 102;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.b0eijwagye
                    section.data(2).logicalSrcIdx = 103;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.oxke0tvvta
                    section.data(3).logicalSrcIdx = 104;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.kr4a3cr0mf
                    section.data(4).logicalSrcIdx = 105;
                    section.data(4).dtTransOffset = 3;

                    ;% rtDW.gdblfyh2ib
                    section.data(5).logicalSrcIdx = 106;
                    section.data(5).dtTransOffset = 4;

                    ;% rtDW.mmq5ttry5y
                    section.data(6).logicalSrcIdx = 107;
                    section.data(6).dtTransOffset = 5;

                    ;% rtDW.agwigalrld
                    section.data(7).logicalSrcIdx = 108;
                    section.data(7).dtTransOffset = 6;

                    ;% rtDW.fudyiowmus
                    section.data(8).logicalSrcIdx = 109;
                    section.data(8).dtTransOffset = 7;

                    ;% rtDW.ftrjzsptcb
                    section.data(9).logicalSrcIdx = 110;
                    section.data(9).dtTransOffset = 8;

                    ;% rtDW.ciu1ngvzgz
                    section.data(10).logicalSrcIdx = 111;
                    section.data(10).dtTransOffset = 9;

                    ;% rtDW.csdsboi20g
                    section.data(11).logicalSrcIdx = 112;
                    section.data(11).dtTransOffset = 10;

                    ;% rtDW.pxsktqyoas
                    section.data(12).logicalSrcIdx = 113;
                    section.data(12).dtTransOffset = 11;

                    ;% rtDW.bya2xiv5j3
                    section.data(13).logicalSrcIdx = 114;
                    section.data(13).dtTransOffset = 12;

                    ;% rtDW.hnl0z3aet1
                    section.data(14).logicalSrcIdx = 115;
                    section.data(14).dtTransOffset = 13;

                    ;% rtDW.ecvugwr44v
                    section.data(15).logicalSrcIdx = 116;
                    section.data(15).dtTransOffset = 14;

                    ;% rtDW.dfdi1srria
                    section.data(16).logicalSrcIdx = 117;
                    section.data(16).dtTransOffset = 15;

                    ;% rtDW.ibwargsoog
                    section.data(17).logicalSrcIdx = 118;
                    section.data(17).dtTransOffset = 16;

                    ;% rtDW.cb3tpkjx5q
                    section.data(18).logicalSrcIdx = 119;
                    section.data(18).dtTransOffset = 17;

                    ;% rtDW.fg4y55crvv
                    section.data(19).logicalSrcIdx = 120;
                    section.data(19).dtTransOffset = 18;

                    ;% rtDW.h4wxpouwz2
                    section.data(20).logicalSrcIdx = 121;
                    section.data(20).dtTransOffset = 19;

                    ;% rtDW.loncelqhha
                    section.data(21).logicalSrcIdx = 122;
                    section.data(21).dtTransOffset = 20;

                    ;% rtDW.afd0kmksnl
                    section.data(22).logicalSrcIdx = 123;
                    section.data(22).dtTransOffset = 21;

                    ;% rtDW.pv0jmcqzdg
                    section.data(23).logicalSrcIdx = 124;
                    section.data(23).dtTransOffset = 22;

                    ;% rtDW.lbzvytb4d4
                    section.data(24).logicalSrcIdx = 125;
                    section.data(24).dtTransOffset = 23;

                    ;% rtDW.pagreemrym
                    section.data(25).logicalSrcIdx = 126;
                    section.data(25).dtTransOffset = 24;

                    ;% rtDW.kkiofr3le1
                    section.data(26).logicalSrcIdx = 127;
                    section.data(26).dtTransOffset = 25;

                    ;% rtDW.pqnrnlrnbt
                    section.data(27).logicalSrcIdx = 128;
                    section.data(27).dtTransOffset = 26;

                    ;% rtDW.d1ahpk1qqr
                    section.data(28).logicalSrcIdx = 129;
                    section.data(28).dtTransOffset = 27;

                    ;% rtDW.hw5n3aggbh
                    section.data(29).logicalSrcIdx = 130;
                    section.data(29).dtTransOffset = 28;

                    ;% rtDW.kgduzxmrs5
                    section.data(30).logicalSrcIdx = 131;
                    section.data(30).dtTransOffset = 29;

            nTotData = nTotData + section.nData;
            dworkMap.sections(5) = section;
            clear section


            ;%
            ;% Non-auto Data (dwork)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        dworkMap.nTotData = nTotData;



    ;%
    ;% Add individual maps to base struct.
    ;%

    targMap.paramMap  = paramMap;
    targMap.signalMap = sigMap;
    targMap.dworkMap  = dworkMap;

    ;%
    ;% Add checksums to base struct.
    ;%


    targMap.checksum0 = 343361676;
    targMap.checksum1 = 2573916448;
    targMap.checksum2 = 664120902;
    targMap.checksum3 = 3813395425;

