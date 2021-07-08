classdef ESDIRK69Const < DI_RKConst
    
    properties
        graphLineStyle = {};
        eval_RHS = true;
    end
    
    properties(SetAccess = protected)
        name  = 'ESDIRK6(5)9[2]SA'
        description = '6th-order, 9-stage, stiffly-accurate, L-Stable ESDIRK method ESDIRK6(5)9[2]SA from C.A. Kennedy, M. H. Carpenter, "Diagonally implicit Runge–Kutta methods for stiff ODEs," 2019';
        order = 5; 
        A = [0                                  0                                   0                                   0                                   0                                   0                                   0                                   0                                   0;
            2/9                                 2/9                                 0                                   0                                   0                                   0                                   0                                   0                                   0;
            1/9                                 -52295652026801/1014133226193379    2/9                                 0                                   0                                   0                                   0                                   0                                   0;
            37633260247889/456511413219805      -162541608159785/642690962402252    186915148640310/408032288622937     2/9                                 0                                   0                                   0                                   0                                   0;
            -37161579357179/532208945751958     -211140841282847/266150973773621    884359688045285/894827558443789     845261567597837/1489150009616527    2/9                                 0                                   0                                   0                                   0;
            32386175866773/281337331200713      498042629717897/1553069719539220    -73718535152787/262520491717733     -147656452213061/931530156064788    -16605385309793/2106054502776008    2/9                                 0                                   0                                   0;
            -38317091100349/1495803980405525    233542892858682/880478953581929     -281992829959331/709729395317651    -52133614094227/895217507304839     -9321507955616/673810579175161      79481371174259/817241804646218      2/9                                 0                                   0;
            -486324380411713/1453057025607868   -1085539098090580/1176943702490991  370161554881539/461122320759884     804017943088158/886363045286999     -15204170533868/934878849212545     -248215443403879/815097869999138    339987959782520/552150039467091     2/9                                 0;
            0                                   0                                   0                                   281246836687281/672805784366875     250674029546725/464056298040646     88917245119922/798581755375683      127306093275639/658941305589808     -319515475352107/658842144391777    2/9
            ];        
       b = [0       0       0                                       281246836687281/672805784366875     250674029546725/464056298040646     88917245119922/798581755375683   127306093275639/658941305589808    -319515475352107/658842144391777    2/9];
       c = [0       4/9     376327483029687/1335600577485745        433625707911282/850513180247701     183/200                             62409086037595/296036819031271   81796628710131/911762868125288     97/100                              1];
    end
    
    methods
        
        function this = ESDIRK69Const(options)
            if(nargin == 0)
                options = struct();
            end
            this@DI_RKConst(options);
        end
    end
    
end