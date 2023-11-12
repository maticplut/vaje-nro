
tocke[n_]:= RandomReal[{-1,1},{n,2}];



preveri[t_]:=Module[{zun={},not={},n=Length[t],i},
i=1;
While[i <=n,
If[t[[i]][[1]]^2+t[[i]][[2]]^2<=1,
not=Append[not,t[[i]]],
zun=Append[zun,t[[i]]]
];
i++;
];
{not,zun}
];


izracunPi[{not_, zun_}] := Module[{priblizekPi, napaka},
   priblizekPi = N[(4 Length[not])/(Length[not] + Length[zun])];
   napaka = Abs[priblizekPi - Pi];
   {priblizekPi, napaka}
   ];

calcPi[] := Module[{i, skPrib = {}, skNpk = {}},
   i = 500;
   While[i <= 20000, Module[{prib, npk},
     {prib, npk} = izracunPi[preveri[tocke[i]]];
     skPrib = Append[skPrib, prib];
     skNpk = Append[skNpk, npk];
     ];
    i += 500;
    ];
   {skPrib, skNpk}
   ];
