% BP Algorithm with GAI test in LLR domain
clear;
clc;

rand('state',114514);
randn('state',114514);

% profile on;
EbN0db      =25:1:35;
SampleNum   = Inf;       % max sample no.aaa
delta       = 0.5; %0.55;          % damping factor
isCorr = 0;
rho = 0.3;  % related channel parameter
max_frame   =1000;
ModType 	= 6;

kbestCount=0;


% Detection Parameter Set
TxAntNum    =5;
RxAntNum    =7;
iterNum     =10;% 7;        % iter no.
N 			= 1; 						% Block length
InfoLen 	= N * ModType * TxAntNum; 	% Information length

switch(ModType)
    case 2
        norm_f=1/sqrt(2);
        QAM=4;
        s=[-1 +1]';
        sym = s*norm_f; % normalized symbol set
        %cons = getConstellation(ModType);
    case 4
        norm_f = 1/sqrt(10);
        QAM=16;
        s=[-1 -3 +1 +3]';
        sym = s*norm_f; % normalized symbol set
        %cons = getConstellation(ModType);
    case 6
        norm_f=1/sqrt(42);
        QAM=64;
        s=[-3 -1 -5 -7 3 1 5 7 ]';
        sym = s*norm_f; % normalized symbol set
        %cons = getConstellation(ModType);
    case 8
        QAM=256;
        norm_f=1/sqrt(170);
        s=[-5 -7 -3 -1 -11 -9 -13 -15 5 7 3 1 11 9 13 15]';
        sym = s*norm_f; % normalized symbol set
        %cons = getConstellation(ModType);
end

slen = length(sym);
H           = zeros(2*RxAntNum,2*TxAntNum);
RxSymbol    = zeros(2*RxAntNum,1); % Received


tic
for nEN =1:length(EbN0db)

    % initialization
    SNR = EbN0db(nEN);
    error_frame = 0;
    error_bits = 0;


    TxBitsLLR = zeros(InfoLen, 1);

    loop = 0;
    while error_frame < max_frame && loop < SampleNum

        loop = loop + 1;

        TxBits = randi([0, 1], [InfoLen, 1]);
        TxSymbol = GrayMapQAM(TxBits', ModType).';
        % TxSymbol=qammod(TxBits,QAM,"gray",'UnitAveragePower',true,'InputType','bit');
        TxSymbol = GrayMapQAM(TxBits', ModType).';
        % TxSymbol=qammod(TxBits,QAM,"gray",'UnitAveragePower',true,'InputType','bit');
        TxSymbol_real=[real(TxSymbol); imag(TxSymbol)];
        % Channel Matrix
        if isCorr
            HMat = channel_mimo(RxAntNum,TxAntNum,rho); %Correlated channel
        else
            HMat = sqrt(0.5) *(randn(RxAntNum, TxAntNum) + sqrt(-1)*randn(RxAntNum, TxAntNum)); %Rayleigh channel
        end




        % Received and sigma^2
        [RxSymbol_c, Nv] = acwgn_EbN0(HMat*TxSymbol, 1, SNR, TxAntNum, RxAntNum, ModType, 1);	% Add Gaussian noise

        if (mod(loop,10) == 0)
            fprintf("\nNow Iter: %d\tNow SNR: %d\tNow FER: %fNow BER: %f\tNow Error Frame: %d \tNow Error Bits: %d", ...
                loop, EbN0db(nEN), error_frame/loop, error_bits/loop/TxAntNum/ModType, error_frame, error_bits);
        end

        % Complex channel matrix to real
        Hreal = [real(HMat), -imag(HMat); imag(HMat), real(HMat)];

        RxSymbol(1:RxAntNum,1)=real(RxSymbol_c);
        RxSymbol(RxAntNum+1:end,1)=imag(RxSymbol_c);

        % symLLR = Bsp_b11_fast(TxAntNum,RxAntNum,slen,RxSymbol,Hreal,Nv,sym,delta,iterNum);[~,Indice] = max(symLLR);symest = sym(Indice);

        % symest = RD_GAI_BP_Det(TxAntNum,RxAntNum,slen,RxSymbol,Hreal,Nv,sym,delta,iterNum);

        % symest = AMP(TxAntNum,RxAntNum,slen,RxSymbol,Hreal,Nv,sym,delta,iterNum);
        % symest= EP(TxAntNum,RxAntNum,slen,RxSymbol,Hreal,Nv,sym,delta,iterNum);
        % symest = GD(TxAntNum,RxAntNum,RxSymbol,Hreal,sym,iterNum,Nv);

        %symest = MDPI_GAI(TxAntNum,RxAntNum,slen,RxSymbol,Hreal,Nv,sym,delta,iterNum);
        symest = test_GAI(TxAntNum,RxAntNum,slen,RxSymbol,Hreal,Nv,sym,delta,iterNum,TxSymbol_real);
        %symest = AltMin(TxAntNum,RxAntNum,slen,RxSymbol,Hreal,Nv,sym,delta,iterNum,TxSymbol_real);
        % symest = MMSE(TxAntNum,RxAntNum,slen,RxSymbol,Hreal,Nv,sym,delta,iterNum);
        % symest =NA(RxSymbol_c,HMat,Nv,sym,norm_f,TxSymbol);
   

        %symest = FC_GAI_BP_Det(TxAntNum,RxAntNum,slen,RxSymbol,Hreal,Nv,sym,delta,iterNum);

        %symest =  Repeat_GAI(TxAntNum,RxAntNum,slen,RxSymbol,Hreal,Nv,sym,delta,iterNum,TxSymbol_real,norm_f);
        %symest =  Repeat_GAI_IDD(TxAntNum,RxAntNum,slen,RxSymbol,Hreal,Nv,sym,delta,iterNum,TxSymbol_real,norm_f);


        % no use
        % symest =  Repeat_Bsp(TxAntNum,RxAntNum,slen,RxSymbol,Hreal,Nv,sym,delta,iterNum,TxSymbol_real,norm_f);
        % symLLR = RD_BsP_df1dm1_Det(TxAntNum,RxAntNum,slen,RxSymbol,Hreal,Nv,sym,delta,iterNum);[~,Indice] = max(symLLR);symest = sym(Indice);

         % DANGER!DANGER!
%         symest=  maxlikehood(RxSymbol,Hreal,sym);


        % QR for K_BEST
        % [bQ, bR] = qr(Hreal);
        % R = bR(1:2*TxAntNum, :);
        % Q = bQ(:, 1:2*TxAntNum);
        % z = Q' * RxSymbol;

%          [symest,kbesttime,PEDtemp] =mcts_det(R, sym, z,Nv,TxSymbol_real); kbestCount=kbestCount+kbesttime;averagekbest=kbestCount/loop;PEDCount=[PEDCount;PEDtemp];
       % symest=K_Best(R, sym', z,128);

        % 用那个祖传格雷码的时候用这4行
                [~,indiceest]=min(abs(sym'-symest),[],2);
                TxBits_est=de2bi(indiceest-1,log2(size(sym,1)),'left-msb');
                TxBits_est=[TxBits_est(1:TxAntNum,:) TxBits_est(TxAntNum+1:end,:)]';
                TxBits_est=TxBits_est(:);
                [~,indiceest]=min(abs(sym'-symest),[],2);
                TxBits_est=de2bi(indiceest-1,log2(size(sym,1)),'left-msb');
                TxBits_est=[TxBits_est(1:TxAntNum,:) TxBits_est(TxAntNum+1:end,:)]';
                TxBits_est=TxBits_est(:);


        % RxSymbol_est=symest(1:TxAntNum)+1j*symest(TxAntNum+1:end);
        % TxBits_est=qamdemod(RxSymbol_est,QAM,"gray","OutputType","bit","UnitAveragePower",true);
        % RxSymbol_est=symest(1:TxAntNum)+1j*symest(TxAntNum+1:end);
        % TxBits_est=qamdemod(RxSymbol_est,QAM,"gray","OutputType","bit","UnitAveragePower",true);

       err=nnz(TxBits_est-TxBits);
       
            % H_error=cat(3,H_error,Hreal);
            % Tx_error=[Tx_error,TxSymbol_real];
            % Rx_error=[Rx_error,RxSymbol];
            % bitsNum_error=[bitsNum_error;err];
            % matlabEst=[matlabEst,(Hreal'*Hreal+Nv*eye(2*TxAntNum)) \ Hreal'*RxSymbol];


        

        
        error_bits = error_bits + err;
        error_frame = error_frame + ~~err;


    end
   
    BER(nEN) = error_bits / TxAntNum / ModType / loop;
    FER(nEN) = error_frame / loop;

end
toc

% writeNPY(H_error,"./npy-matlab/H_error.npy")
% writeNPY(Tx_error,"./npy-matlab/Tx_error.npy")
% writeNPY(Rx_error,"./npy-matlab/Rx_error.npy")
% writeNPY(bitsNum_error,"./npy-matlab/bitsNum_error.npy")
% writeNPY(matlabEst,"./npy-matlab/matlabEst.npy")
% writeNPY(H_error,"./npy-matlab/H_error.npy")
% writeNPY(Tx_error,"./npy-matlab/Tx_error.npy")
% writeNPY(Rx_error,"./npy-matlab/Rx_error.npy")
% writeNPY(bitsNum_error,"./npy-matlab/bitsNum_error.npy")
% writeNPY(matlabEst,"./npy-matlab/matlabEst.npy")
