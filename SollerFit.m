%根据目标文件中的荧光信息，拟合出具体x，y，z的位置
%输入参数：filepath、zposi
%输出参数：xfit、yfit、zfit、vari

function [xfit, yfit, zfit, xPosi, isfit, origSollervalue] = SollerFit(filepath, zposi, chnlmin, chnlmax)
    %初始参数：N->狭缝数量，B->下底长度，D->上底长度,H->狭缝高度，L->狭缝焦距
    N = 7;
    B = 57.417;
    D = 28.7;
    H = 50;
    L = 100;
    %设置拟合收敛条件
    
    %读取filepath中的Soller狭缝信息
    [Soller_Value, xPosi] = GetData(filepath, chnlmin, chnlmax);
    origSollervalue = Soller_Value;
    %disp(Soller_Value);
    %消去奇异点
    Soller_Value = delerrorpoint(Soller_Value, N - 1);
    
    %判断狭缝是否可以被拟合
    [whethcbfitflag, movedir] = whethcbfit(Soller_Value, N - 1);
    if whethcbfitflag == 1
        disp('================================================================');
        %拟合狭缝结果
        [ result ] = nfitsoller(Soller_Value, 1, 0);
        %disp('result');disp(result(2:3));
        %求解z的值
        zfit = calczposi(result(1), result(4), N, H, B, D);
        %求解x，y的值
        xfit = calcxposi(result(1), result(2), zfit, N, H, B, D);
        yfit = calcxposi(result(1), result(3), zfit, N, H, B, D);
        isfit = 1;
    elseif whethcbfitflag == 0
        disp('================================================================');
        disp('数据不可拟合，需提前移动')
        %拟合狭缝结果
        nfitsoller(Soller_Value, 1, 0);
        [result] = estimateposition(Soller_Value, zposi, N, H, B, D, 1);
        %disp('result');disp(result);
        %求解x，y的值
        k = 2*((B-D)*zposi - B*H)/(2*(N-1)*(zposi-H));
        xfit = calcxposi(k*result(1), result(2), zposi, N, H, B, D);
        yfit = calcxposi(k*result(1), result(3), zposi, N, H, B, D);
        zfit = zposi;
        isfit = 0;
    elseif whethcbfitflag == 2
        disp('================================================================');
        disp('数据平均信号小于500，不可被拟合');
        %拟合狭缝结果
        nfitsoller(Soller_Value, 1, 0);
        %求解信号分布位置
        xfit = 0%3 * movedir;
        yfit = 0;
        zfit = zposi;
        isfit = 0;
    elseif whethcbfitflag == 3
        disp('================================================================');
        disp('数据通道小于6，不可被拟合'); 
        %求解信号分布位置
        xfit = 0%2 * movedir;
        yfit = 0;
        zfit = zposi;
        isfit = 0;
    end
    %disp('timestamp:',[]);
    disp('Soller fit result is:');
    disp(['x:',num2str(xfit)]);
    disp(['y:',num2str(yfit)]);
    disp(['z:',num2str(zfit)]);
    disp('================================================================');
    disp(char(10));
end

%======================================================================================%
function S = delerrorpoint(S, N)
    %disp(S);
    S2 = [S(:,:)];
    stdS2 = std(S2(S2 ~= 0));
    meanS2 = mean(S2(S2 ~=0));
    %disp(stdS2);
    S2=S2 .* (abs(S2-meanS2) <= 4*stdS2);
    S = reshape(S2, N, N);
    %disp(S);
end


%======================================================================================%
function [rslt dir] = whethcbfit(S, N)
    rslt = 0;
    %判断荧光信号平均数
    Scut = S(:);
    Scut = Scut(Scut ~= 0);
    %如果每个通道荧光信号平均数小于500，则不可被拟合
    if mean(Scut) < 500
        %disp('singal number less than 500 per channel!');
        rslt = 2;
    end
    %如果有计数的荧光通道小于6，则不可被拟合
    if numel(Scut) < 6
        %disp('fitable channel are less than 6!');
        rslt = 3;
    end
    %找到y的最大值
    [maxs,mp] = max(S);
    maxyposi = mode(mp);
    mpy = mp;
    %S进行转置，求解x的最大值
    S = S.';
    [maxs,mp] = max(S);
    maxxposi = mode(mp);
    mpx = mp;
    
%     disp('maxxposi');
%     disp(maxxposi);
%     disp('maxyposi');
%     disp(maxyposi);
    if rslt == 0
        if maxyposi == N || maxxposi == N || maxyposi == 1 || maxxposi == 1
            rslt = 0;
        else
            rslt = 1;
        end
    end
    dir = sign(maxxposi - 3);
end

%======================================================================================%
function [ result ] = nfitsoller(S, elimzero, dispfig)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
    %S是6*6个通道的荧光数据
    %行为x,列为y
    %elimzero 为是否剔除零点
    %dispfig 是否显示图片
%     S=[ 0       2874	5317	0       0       0
%         0       7050	15139	19316	19452	11397
%         6302	16964	28947	35778	34690	12569
%         8170	23596	0       39304	47513	27921
%         6386	19919	34227	46017	40252	29449
%         0       10009	19043	0       0       0];
% S=[        4608        7744        9472        7808        5120        2368
%         7992       13431       16428       13542        8880        4107
%        10512       17666       21608       17812       11680        5402
%         9648       16214       19832       16348       10720        4958
%         6552       11011       13468       11102        7280        3367
%         3240        5445        6660        5490        3600        1665];

    %找到y的最大值作
    [maxs,mp] = max(S);
    maxyposi = mode(mp);
    mpy = mp;
    %S进行转置，求解x的最大值
    S = S.';
    [maxs,mp] = max(S);
    maxxposi = mode(mp);
    mpx = mp;
    %找到S的最大值作为直线最高点初值
    maxs = max(maxs);
    if dispfig == 1
        %在一副图里面绘制S的分布
        figure(1);
        bar3(S);
    end
	
    %构造需要拟合的数组
    %x轴123..123..123....
    x = repmat([1 : 1 : 6], 1, 6);
    %y轴111..222..333....
    y = repmat([1 : 1 : 6], 6, 1);
    y = y(:)';
    %确定k的符号, max点往左为上升趋势，max往右为下降趋势
    %max的符号无影响仅是确定初始值使用
    xf = sign(maxxposi - x + 0.5);
    yf = sign(maxyposi - y + 0.5);
    sourceS = S(:)';
    %拟合函数的四个已知量x，y，xf，yf和未知量S
    if elimzero == 1
        %去除掉零点和无效点
        x(find(sourceS==0))=[];
        y(find(sourceS==0))=[];
        xf(find(sourceS==0))=[];
        yf(find(sourceS==0))=[];
        sourceS(find(sourceS==0))=[];
    end
    sourcex = [x' y' xf' yf'];
    sourcex = sourcex';
    
    %f=@(a,x)(a(1)*x(3,:)*(x(1,:) - a(2)) + a(4))*(a(1)*x(4,:)*(x(2,:) - a(3)) + a(4));
    %----k = a1, a1 = a2, a2 = a3, b = a4; 
    %确定初始值
    %按照y方向进行拟合
    for i = 1 : 1 : 6
        if(maxxposi >= 3)
            lx = 1 : 1 : mpy(i);
            ly = S(i, 1 : mpy(i));
        else
            lx = mpy(i) : 1 : 6;
            ly = S(i, mpy(i) : 6);
        end
        rf = polyfit(lx, ly, 1); 
        rk(i) = rf(1);
    end
    
    a(1) = sign(rk(maxxposi))*(rk(maxxposi)^2)^0.25;
    a(2)= maxxposi;
    a(3) = maxyposi;
    a(4) = maxs;
    
    initvalue = a;
    
    option=optimset('MaxFunEvals',2^12,'MaxIter',2^14,'TolX',1e-10,'TolFun',1e-10);
    
    [result,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@fun, initvalue, sourcex, sourceS, [], [], option);
end

%======================================================================================%
function [ f ] = fun( a, x )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%f=1-8/9.8696.*exp(-9.8696.*c(1).*x/(c(2).^2));
    f=(a(1)*x(3,:).*(x(1,:) - a(2)) + a(4)).*(a(1)*x(4,:).*(x(2,:) - a(3)) + a(4));
end

%======================================================================================%
function [Soller_Value, xPos] = GetData(filepath, chnlmin, chnlmax)
    dispflag = 0;
    if(isempty(dispflag))
        filepath = 'E:\工作文档\XAFS\课题\soller slits\Model\realdata\Fe_70\Fe_Z70_XN32_R1';
        chnlmin = 1;
        chnlmax = 0;
        dispflag = 1;
    end
    fid = fopen(filepath);

    
    while ~feof(fid)
        str = fgetl(fid);
        if(strfind(str, 'X position:') > 0)
            buff = strsplit(str,':');
            buff = char(buff(2));
            %disp(char(buff))
            xPos = buff;
        end
        if(strfind(str, 'PV Number:') > 0)
            data.PV_Number = '';
            data.Channel = '';
            data.R0LO = '';
            data.R0HI = '';
            data.R0 = '';
            data.MCABinWidth_RBV = '';
            data.Timestamp = '';
            data.Data = '';
            buff = strsplit(str,':');
            buff = char(buff(2));
            data.PV_Number = buff;
            while(isempty(strfind(str,'=============================')) == 1)
                str = fgetl(fid);
                %disp(str);
                if(strfind(str, 'Channel:') > 0)
                    buff = strsplit(str,':');
                    buff = char(buff(2));
                    %disp(char(buff))
                    data.Channel = buff;
                    if(str2num(buff) == 36)
                        %disp('OK!');
                        break;
                    end
                elseif(strfind(str, 'R0LO:') > 0)
                    buff = strsplit(str,':');
                    buff = char(buff(2));
                    data.R0LO = buff;
                elseif(strfind(str, 'R0HI:') > 0)
                    buff = strsplit(str,':');
                    buff = char(buff(2));
                    data.R0HI = buff;
                elseif(strfind(str, 'R0:') > 0)
                    buff = strsplit(str,':');
                    buff = char(buff(2));
                    data.R0 = buff;
                elseif(strfind(str, 'MCABinWidth_RBV:') > 0)
                    buff = strsplit(str,':');
                    buff = char(buff(2));
                    data.MCABinWidth_RBV = buff;
                elseif(strfind(str, 'Timestamp:') > 0)
                    buff = strsplit(str,':');
                    buff = char(buff(2));
                    data.Timestamp = buff;
                elseif(strfind(str, 'Data:') > 0)
                    str = fgetl(fid);
                    data.Data = char(str);
    %                 disp(data.Data);
                end
            end
            id = int8(str2num(data.Channel));
            %disp(class(id));
            rslt(id) = data;
        end
    end
    
    if chnlmin ~= 0 || chnlmax ~= 0
        chnlmin = chnlmin + 1;
        chnlmax = chnlmax + 1;
        for i = 1 : 1 : 6
            for j = 1 : 1 : 6
                data = rslt((i-1)*6 + j);
                datastr = data.Data;
                if strcmp(datastr, '') == 0
                    expr = '[^\t]*[^\t]*';
                    alldata = regexp(datastr, expr,'match');
                    alldata = char(alldata);
                    alldata = str2num(alldata);
                    Soller_Value(i , j) = sum(alldata(chnlmin:chnlmax));
                else
                    Soller_Value(i , j) = 0; 
                end 
            end
        end
    else
        for i = 1 : 1 : 6
            for j = 1 : 1 : 6
                data = rslt((i-1)*6 + j);
        %         disp(str2num(data.R0))
        %         disp(size(str2num(data.R0)))
                Soller_Value(i , j) = str2double(data.R0);
                %disp(map(2, (i-1)*6 + j));
            end
        end
    end

    %disp(size(rslt))
    if(dispflag)

        disp(Soller_Value);
        bar3(Soller_Value);
        box off;

        for j = 1 : 1 : 6
            for i = 1 : 1 : 6
                text(j,i,round(Soller_Value(i,j)), num2str(round(Soller_Value(i,j))),'fontsize',10)
            end
        end
        xlabel('class','fontsize',12,'fontname','times')
        ylabel('record','fontsize',12,'fontname','times')
        zlabel('frequence','fontsize',12,'fontname','times')
    end
    fclose(fid);

    Soller_Value(find(isnan(Soller_Value) == 1)) = 0;
end

%======================================================================================%
function [z, q] = calczposi(k, b, N, H, B, D)
    z = (-2*B*H*b + k*B*H)/(k*(B+D) - 2*b*(B-D));
%     func1 = [num2str(k*2*(N-1)), '*q*(z-', num2str(H), ')-2*(',num2str(B-D),'*z+',num2str(B*H),'=0']
%     func2 = [num2str(b*2*(N-1)), '*q*(z-', num2str(H), ')-',num2str(D + B),'*z-',num2str(B*H),'=0']
%    if num2str(sign(2*((B-D)*z-B*H))) == 1
%        sig = '+';
%    else
%        sig = '-';
%    end
%    func1 = [num2str(k*2*(N-1)*(z-H)), '*q', sig, num2str(abs(2*((B-D)*z-B*H))),'=0'];
%    q = solve(func1);
end

%======================================================================================%
function x = calcxposi(k, a, z, N, H, B, D)
%     x = ((N+1)*(z-H)*B + B*H - k*a - (N+1)*D*z - (B - D)*z)/(2*(N-1)*H);
      %disp(k);
      k = 2*((B-D)*z - B*H);
      %disp(k);
      x = ((B-D)*N*z - B*N*H - k*a)/(2*(N-1)*H);
end

%======================================================================================%
function [result] = estimateposition(Soller_Value, z, N, H, B, D, elimzero)
    k = 2*((B-D)*z - B*H);%/(2*(N-1)*(z-H));
    b = ((B+D)*z - B*H);%/(2*(N-1)*(z-H));
    S = Soller_Value;
%     disp('k');   disp(k);
%     disp('b');   disp(b);

    %找到y的最大值作
    [maxs,mp] = max(S);
    maxyposi = mode(mp);
    mpy = mp;
    %S进行转置，求解x的最大值
    S = S.';
    [maxs,mp] = max(S);
    maxxposi = mode(mp);
    mpx = mp;
    %找到S的最大值作为直线最高点初值
    maxs = max(maxs);
    
    %构造需要拟合的数组
    %x轴123..123..123....
    x = repmat([1 : 1 : 6], 1, 6);
    %y轴111..222..333....
    y = repmat([1 : 1 : 6], 6, 1);
    y = y(:)';
    %确定k的符号, max点往左为上升趋势，max往右为下降趋势
    %max的符号无影响仅是确定初始值使用
    xf = sign(maxxposi - x + 0.5);
    yf = sign(maxyposi - y + 0.5);
    sourceS = S(:)';
    %拟合函数的四个已知量x，y，xf，yf和未知量S
    if elimzero == 1
        %去除掉零点和无效点
        x(find(sourceS==0))=[];
        y(find(sourceS==0))=[];
        xf(find(sourceS==0))=[];
        yf(find(sourceS==0))=[];
        sourceS(find(sourceS==0))=[];
    end
    sourcex = [x' y' xf' yf'];
    sourcex = sourcex';    
    %确定初始值
    %按照y方向进行拟合
    for i = 1 : 1 : 6
        if(maxxposi >= 3)
            lx = 1 : 1 : mpy(i);
            ly = S(i, 1 : mpy(i));
        else
            lx = mpy(i) : 1 : 6;
            ly = S(i, mpy(i) : 6);
        end
        rf = polyfit(lx, ly, 1); 
        rk(i) = rf(1);
    end
    
    
    a(1) = mean(rk)/k;
    a(2) = maxxposi;
    a(3) = maxyposi;
    
    initvalue = a;
    
    fun2 = @(a, x)(a(1)^2*(k * x(3, :) .* (x(1,:) - a(2)) + b) .* (k * x(4, :) .* (x(2,:) - a(3)) + b));
    
    option=optimset('MaxFunEvals',2^12,'MaxIter',2^14,'TolX',1e-10,'TolFun',1e-10);
    
    [result,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(fun2, initvalue, sourcex, sourceS, [], [], option);
end
