function piv_frame_browser(ii)
% piv_frame_browser  PIVlab-style frame browser for longitudinal PIV.
%
%   piv_frame_browser(ii)  — experiment index ii (default 4).
%
%   Left panel : raw PIV image.  Flip A/B [space], draggable IW box,
%                displacement calculator (pixel coords A/B).
%   Right panel: velocity field, populated ONLY when you click "Show ->".
%                Changing the frame or source does NOT auto-update the right
%                panel — click Show -> to compute/load for the current frame.
%
%   Velocity source (bottom row):
%     "Precomputed file"    : load Chris_recompute/PIVMat/*_compVel_*.mat
%     "Recompute (rect IW)" : run ComputeVelocities on the editable IW_x/IW_z
%                              pyramids with 50% or 75% overlap.
%                              Nothing saved to disk.
%
%   Navigation: slider, Prev/Next, frame edit box, left/right arrow keys.
%   Image data cached per frame.  Velocity cached per frame+pyramid+overlap.

if nargin < 1 || isempty(ii), ii = 4; end

%% ------------------------------------------------------------------------
%% Config
%% ------------------------------------------------------------------------
cfg.LONG     = 'D:\DelawareDataBackup\Longitudinal\PIV\';
cfg.rootpath = 'D:\Scripps';
cfg.DX = 1/17697.69;
cfg.DT = 10e-3;
cfg.IntrWndw      = [256 128 64 32 16 8];   % assumed square pyramid for precomputed display
cfg.num_of_digits = 3;
cfg.clim_u = [-0.01, 0.12];
cfg.clim_w = [-0.04, 0.04];
cfg.fieldLabels = {'dx raw (px/fr)','dx smooth (px/fr)', ...
                   'dz raw (px/fr)','dz smooth (px/fr)', ...
                   'u (m/s)','w (m/s)','dcor'};
cfg.fieldKeys   = {'dx_raw','dx_smooth','dz_raw','dz_smooth','u','w','dcor'};

addpath(strcat(cfg.rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\'));
addpath(strcat(cfg.rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\FabriceScripts\'));
addpath(strcat(cfg.rootpath,'\GC-Wave-Gen\M-Files_FabMarcNovDec2014\CrapperOptimizedFindSurface\'));

DIRS = dir(cfg.LONG); DIRS = DIRS(3:end);
assert(ii >= 1 && ii <= numel(DIRS), 'ii=%d out of range (1..%d)', ii, numel(DIRS));
cfg.exp_name  = DIRS(ii).name;
cfg.load_path = [cfg.LONG cfg.exp_name];

fa   = dir([cfg.load_path filesep 'PIVRaw' filesep 'PIV' filesep ...
            cfg.exp_name '_Piv_*_a.mat']);
tok  = regexp({fa.name}, '_Piv_(\d+)_a\.mat$', 'tokens', 'once');
keep = ~cellfun('isempty', tok);
frames = sort(cellfun(@(c) str2double(c{1}), tok(keep)));
assert(~isempty(frames), 'No PIV frames found under %s', cfg.load_path);

try
    cmap = brewermap([], 'Spectral');
catch
    cmap = parula(256);
end

%% ------------------------------------------------------------------------
%% Shared state
%% ------------------------------------------------------------------------
S             = struct();
S.cfg         = cfg;
S.frames      = frames;
S.idx         = find(frames == 157, 1);
if isempty(S.idx), S.idx = 1; end
S.ab          = 'a';
S.field       = 'dx_raw';
S.quiver      = false;
S.qscale      = 1;
S.recompute   = false;              % false = load precomputed file
S.IntrWndw_x  = [64 32 16 8];      % recompute x-pyramid (IW width)
S.IntrWndw_z  = [16  8  8 4];      % recompute z-pyramid (IW height, rect)
S.overlapFrac = 0.5;               % 0.50 = 50%  |  0.25 = 75% overlap
S.GrdSpc_x    = S.IntrWndw_x / 2;
S.GrdSpc_z    = S.IntrWndw_z / 2;
% IW box size — starts from precomputed square pyramid (display only)
S.iw_x        = cfg.IntrWndw(end);
S.iw_z        = cfg.IntrWndw(end);
S.roiPos      = [];
S.cmap        = cmap;
S.imageCache  = containers.Map('KeyType','char','ValueType','any');  % images only
S.rawCache    = containers.Map('KeyType','char','ValueType','any');  % velocity
S.imgs        = [];    % current-frame images (left panel)
S.data        = [];    % last computed velocity (right panel, may lag S.idx)
S.ptA         = [NaN NaN];
S.ptB         = [NaN NaN];
S.hImgL  = []; S.hSurfL = []; S.roi = [];
S.hPtA   = []; S.hPtB   = [];
S.hImgR  = []; S.hQuiver = [];

%% ------------------------------------------------------------------------
%% Build UI   (figure 1440 x 880, five control rows)
%% ------------------------------------------------------------------------
S.fig = figure('Name','PIV frame browser','NumberTitle','off', ...
    'Color','white','Position',[50 50 1440 880], ...
    'KeyPressFcn',@onKey);

S.axL = axes('Parent',S.fig,'Position',[0.04 0.31 0.43 0.65]);
S.axR = axes('Parent',S.fig,'Position',[0.54 0.31 0.43 0.65]);

ctl = @(varargin) uicontrol('Parent',S.fig,'Units','normalized',varargin{:});

%% Row 1 (y~0.250): frame label + go-to edit
S.uiFrameLbl = ctl('Style','text','String','', ...
    'Position',[0.04 0.250 0.28 0.033], ...
    'HorizontalAlignment','left','BackgroundColor','white','FontWeight','bold');
ctl('Style','text','String','Go to frame:', ...
    'Position',[0.33 0.250 0.09 0.033], ...
    'HorizontalAlignment','right','BackgroundColor','white');
S.uiFrameEdit = ctl('Style','edit','String','', ...
    'Position',[0.43 0.252 0.06 0.033],'Callback',@onFrameEdit);

%% Row 2 (y~0.192): Prev / slider / Next
N = numel(frames);
S.uiPrev = ctl('Style','pushbutton','String','< Prev', ...
    'Position',[0.04 0.190 0.06 0.045],'Callback',@(~,~) stepFrame(-1));
S.uiNext = ctl('Style','pushbutton','String','Next >', ...
    'Position',[0.43 0.190 0.06 0.045],'Callback',@(~,~) stepFrame(+1));
S.uiSlider = ctl('Style','slider','Min',1,'Max',max(N,2),'Value',S.idx, ...
    'SliderStep',[1/max(N-1,1) min(10/max(N-1,1),1)], ...
    'Position',[0.11 0.190 0.31 0.045],'Callback',@(h,~) setIdx(round(h.Value)));

%% Row 3 (y~0.130): Flip A/B + IW popup  |  Field popup + Show button
S.uiBlink = ctl('Style','pushbutton','String','Flip A/B  [space]', ...
    'Position',[0.04 0.130 0.13 0.045],'FontWeight','bold','Callback',@(~,~) toggleAB());
S.uiAbLbl = ctl('Style','text','String','Viewing: A', ...
    'Position',[0.18 0.130 0.09 0.045], ...
    'FontWeight','bold','ForegroundColor',[0 0.6 0],'BackgroundColor','white');
ctl('Style','text','String','IW (px):', ...
    'Position',[0.27 0.133 0.06 0.036], ...
    'HorizontalAlignment','right','BackgroundColor','white');
S.uiIW = ctl('Style','popupmenu','String',iwStrList(), ...
    'Value',numel(dispIW_x()), ...
    'Position',[0.34 0.135 0.09 0.036],'Callback',@onIW);

ctl('Style','text','String','Field:', ...
    'Position',[0.54 0.133 0.05 0.036], ...
    'HorizontalAlignment','right','BackgroundColor','white');
S.uiField = ctl('Style','popupmenu','String',cfg.fieldLabels,'Value',1, ...
    'Position',[0.60 0.135 0.13 0.036],'Callback',@onField);
S.uiShow = ctl('Style','pushbutton','String','Show  ->', ...
    'Position',[0.74 0.130 0.13 0.045],'FontWeight','bold', ...
    'BackgroundColor',[0.75 1.0 0.75],'Callback',@onShow);

%% Row 4 (y~0.074): displacement calc  |  quiver + scale
ctl('Style','text','String','A  x:', ...
    'Position',[0.04 0.074 0.04 0.036], ...
    'HorizontalAlignment','right','BackgroundColor','white');
S.uiXa = ctl('Style','edit','String','', ...
    'Position',[0.08 0.076 0.045 0.036],'Callback',@onDisp);
ctl('Style','text','String','y:', ...
    'Position',[0.126 0.074 0.02 0.036], ...
    'HorizontalAlignment','right','BackgroundColor','white');
S.uiYa = ctl('Style','edit','String','', ...
    'Position',[0.147 0.076 0.045 0.036],'Callback',@onDisp);
ctl('Style','text','String','B  x:', ...
    'Position',[0.205 0.074 0.04 0.036], ...
    'HorizontalAlignment','right','BackgroundColor','white');
S.uiXb = ctl('Style','edit','String','', ...
    'Position',[0.245 0.076 0.045 0.036],'Callback',@onDisp);
ctl('Style','text','String','y:', ...
    'Position',[0.291 0.074 0.02 0.036], ...
    'HorizontalAlignment','right','BackgroundColor','white');
S.uiYb = ctl('Style','edit','String','', ...
    'Position',[0.312 0.076 0.045 0.036],'Callback',@onDisp);
ctl('Style','text','String','|D| =', ...
    'Position',[0.365 0.074 0.04 0.036], ...
    'HorizontalAlignment','right','BackgroundColor','white');
S.uiDispLbl = ctl('Style','text','String','--', ...
    'Position',[0.405 0.074 0.12 0.036], ...
    'HorizontalAlignment','left','FontWeight','bold','BackgroundColor','white');

S.uiQuiver = ctl('Style','checkbox','String','Quiver', ...
    'Value',0,'Position',[0.54 0.076 0.08 0.036], ...
    'BackgroundColor','white','Callback',@onQuiver);
ctl('Style','text','String','Scale:', ...
    'Position',[0.63 0.076 0.05 0.036], ...
    'HorizontalAlignment','right','BackgroundColor','white');
S.uiQscale = ctl('Style','slider','Min',0.1,'Max',10,'Value',1, ...
    'SliderStep',[0.1/9.9 1/9.9], ...
    'Position',[0.69 0.076 0.20 0.036],'Callback',@onQscale);
S.uiQscaleLbl = ctl('Style','text','String','1.0', ...
    'Position',[0.90 0.076 0.05 0.036], ...
    'HorizontalAlignment','left','BackgroundColor','white');

%% Row 5 (y~0.016): source | IW_x  IW_z  overlap  GrdSpc label
ctl('Style','text','String','Velocity:', ...
    'Position',[0.04 0.014 0.06 0.036], ...
    'HorizontalAlignment','right','BackgroundColor','white');
S.uiSource = ctl('Style','popupmenu', ...
    'String',{'Precomputed file','Recompute (rect IW)'},'Value',1, ...
    'Position',[0.105 0.016 0.135 0.036],'Callback',@onSource);
ctl('Style','text','String','IW x:', ...
    'Position',[0.245 0.014 0.04 0.036], ...
    'HorizontalAlignment','right','BackgroundColor','white');
S.uiIntrWndw_x = ctl('Style','edit','String',num2str(S.IntrWndw_x), ...
    'Enable','off','Position',[0.288 0.016 0.12 0.036],'Callback',@onIntrWndw_x, ...
    'TooltipString','Even non-increasing (e.g. 64 32 16 8)');
ctl('Style','text','String','z:', ...
    'Position',[0.412 0.014 0.025 0.036], ...
    'HorizontalAlignment','right','BackgroundColor','white');
S.uiIntrWndw_z = ctl('Style','edit','String',num2str(S.IntrWndw_z), ...
    'Enable','off','Position',[0.440 0.016 0.09 0.036],'Callback',@onIntrWndw_z, ...
    'TooltipString','Even non-increasing, same # levels as IW_x (e.g. 16 8 8 4)');
S.uiOverlap = ctl('Style','popupmenu','String',{'50% overlap','75% overlap'},'Value',1, ...
    'Enable','off','Position',[0.535 0.016 0.09 0.036],'Callback',@onOverlap);
S.uiGsLbl = ctl('Style','text','String',gsLabel(), ...
    'Position',[0.630 0.014 0.33 0.036], ...
    'HorizontalAlignment','left','BackgroundColor','white');

guidata(S.fig, S);
loadImgAndShowLeft();

%% ========================================================================
%% Nested functions
%% ========================================================================
    function pair = pairStr(fnum)
        pair = sprintf(['%0' num2str(cfg.num_of_digits) 'd'], fnum);
    end

    % Active display IW vectors (precomputed = square cfg pyramid; recompute = rect)
    function v = dispIW_x()
        if S.recompute, v = S.IntrWndw_x; else, v = cfg.IntrWndw; end
    end
    function v = dispIW_z()
        if S.recompute, v = S.IntrWndw_z; else, v = cfg.IntrWndw; end
    end

    function s = iwStrList()
        vx = dispIW_x();  vz = dispIW_z();
        if isequal(vx, vz)
            s = arrayfun(@(v) sprintf('%d',v), vx, 'UniformOutput', false);
        else
            s = arrayfun(@(x,z) sprintf('%dx%d',x,z), vx, vz, 'UniformOutput', false);
        end
    end

    function lbl = gsLabel()
        lbl = sprintf('GrdSpc  x=%s   z=%s', ...
            num2str(S.GrdSpc_x), num2str(S.GrdSpc_z));
    end

    function k = cacheKey(fnum)
        if S.recompute
            k = sprintf('R:%d:x%s_z%s_ov%d', fnum, ...
                num2str(S.IntrWndw_x), num2str(S.IntrWndw_z), ...
                round(S.overlapFrac * 100));
        else
            k = sprintf('P:%d', fnum);
        end
    end

    % ------------------------------------------------------------------
    % loadImages  – load raw images + surface masks; cached per frame only.
    % ------------------------------------------------------------------
    function imgs = loadImages(fnum)
        key = sprintf('IMG:%d', fnum);
        if isKey(S.imageCache, key), imgs = S.imageCache(key); return; end
        pair = pairStr(fnum);
        lp   = cfg.load_path;
        Aa = load([lp filesep 'PIVRaw' filesep 'PIV' filesep ...
                   cfg.exp_name '_Piv_' pair '_a.mat']);
        Bb = load([lp filesep 'PIVRaw' filesep 'PIV' filesep ...
                   cfg.exp_name '_Piv_' pair '_b.mat']);
        imgs.IM_a = Aa.imgPiv;
        imgs.IM_b = Bb.imgPiv;
        sa = FindSurfaceCapillary( ...
            [lp filesep 'PIVRaw' filesep 'PIVSURF' filesep ...
             cfg.exp_name '_Pivsurf_' pair '_a.mat'], findMask=true);
        sb = FindSurfaceCapillary( ...
            [lp filesep 'PIVRaw' filesep 'PIVSURF' filesep ...
             cfg.exp_name '_Pivsurf_' pair '_b.mat'], findMask=true);
        imgs.surf_a = sa.surfacePIVImg;
        imgs.surf_b = sb.surfacePIVImg;
        imgs.mask_a = sa.mask;
        imgs.mask_b = sb.mask;
        S.imageCache(key) = imgs;
    end

    % ------------------------------------------------------------------
    % loadRaw  – load/compute velocity; cached per frame+source+pyramid+overlap.
    %            Calls loadImages internally (uses imageCache if warm).
    % ------------------------------------------------------------------
    function raw = loadRaw(fnum)
        key = cacheKey(fnum);
        if isKey(S.rawCache, key), raw = S.rawCache(key); return; end
        imgs = loadImages(fnum);
        pair = pairStr(fnum);
        lp   = cfg.load_path;
        if S.recompute
            IW = [S.IntrWndw_x(:), S.IntrWndw_z(:)];  % Nx2 rect pyramid
            GS = [S.GrdSpc_x(:),   S.GrdSpc_z(:)];
            compVel = ComputeVelocities_Quick_Filt_Deform_Water_dcorFilt( ...
                imgs.IM_a, imgs.IM_b, imgs.mask_a, imgs.mask_b, IW, GS);
        else
            L = load([lp filesep 'Chris_recompute' filesep 'PIVMat' filesep ...
                      cfg.exp_name '_compVel_' pair '.mat']);
            compVel = L.compVel;
        end
        raw.dx0  = compVel.delta_x .* compVel.Mask;
        raw.dz0  = compVel.delta_z .* compVel.Mask;
        raw.dcor = compVel.dcor;
        raw.Mask = compVel.Mask;
        raw.x    = compVel.xPIV * cfg.DX;   % m
        raw.z    = compVel.zPIV * cfg.DX;   % m
        raw.fnum = fnum;
        S.rawCache(key) = raw;
    end

    % ------------------------------------------------------------------
    % deriveFields  – removeOutliers + smoothn.
    % ------------------------------------------------------------------
    function d = deriveFields(raw)
        dx_r = removeOutliers(raw.dx0, raw.dcor);
        dz_r = removeOutliers(raw.dz0, raw.dcor);
        dx_s = smoothn(dx_r, 0.4);
        dz_s = smoothn(dz_r, 0.4);
        d           = raw;
        d.dx_raw    = raw.dx0;
        d.dz_raw    = raw.dz0;
        d.dx_smooth = dx_s .* raw.Mask;
        d.dz_smooth = dz_s .* raw.Mask;
        d.u         = d.dx_smooth * cfg.DX / cfg.DT;
        d.w         = d.dz_smooth * cfg.DX / cfg.DT;
    end

    % ------------------------------------------------------------------
    % loadImgAndShowLeft  – navigation: images only, no velocity.
    % ------------------------------------------------------------------
    function loadImgAndShowLeft()
        fnum = S.frames(S.idx);
        set(S.fig,'Pointer','watch');
        if ~isKey(S.imageCache, sprintf('IMG:%d',fnum))
            set(S.uiFrameLbl,'String', ...
                sprintf('Loading frame %s ...', pairStr(fnum)));
            drawnow;
        end
        S.imgs = loadImages(fnum);
        set(S.fig,'Pointer','arrow');
        renderLeftFull();
        updateFrameLabel();
        guidata(S.fig, S);
    end

    % ---- left panel ----
    function renderLeftFull()
        cla(S.axL, 'reset');
        if isempty(S.imgs), return; end
        IM = pickRaw();
        S.hImgL  = imagesc(S.axL, IM, [0 300]);
        colormap(S.axL, 'gray');
        hold(S.axL, 'on');
        S.hSurfL = plot(S.axL, pickSurf(), '-r', 'LineWidth', 1);

        ptAx = NaN; ptAy = NaN; ptBx = NaN; ptBy = NaN;
        if ~any(isnan(S.ptA)), ptAx = S.ptA(1); ptAy = S.ptA(2); end
        if ~any(isnan(S.ptB)), ptBx = S.ptB(1); ptBy = S.ptB(2); end
        S.hPtA = plot(S.axL, ptAx, ptAy, '+c', 'MarkerSize',14,'LineWidth',2);
        S.hPtB = plot(S.axL, ptBx, ptBy, '+m', 'MarkerSize',14,'LineWidth',2);

        axis(S.axL,'tight'); axis(S.axL,'equal');
        title(S.axL, leftTitle(), 'Interpreter','none');

        sz = size(IM);
        if isempty(S.roiPos)
            pos = [sz(2)/2 - S.iw_x/2, sz(1)/2 - S.iw_z/2, S.iw_x, S.iw_z];
        else
            pos = [S.roiPos(1), S.roiPos(2), S.iw_x, S.iw_z];
        end
        S.roi = drawrectangle(S.axL, 'Position', pos, 'Color', 'y', ...
            'LineWidth', 1, 'InteractionsAllowed', 'translate', ...
            'Label', sprintf('IW %dx%d px', S.iw_x, S.iw_z), 'FaceAlpha', 0);
        addlistener(S.roi, 'ROIMoved', @(src,~) onROIMoved(src));
        hold(S.axL, 'off');
    end

    function renderLeftBlink()
        if isempty(S.hImgL) || ~isvalid(S.hImgL), renderLeftFull(); return; end
        S.hImgL.CData  = pickRaw();
        S.hSurfL.YData = pickSurf();
        title(S.axL, leftTitle(), 'Interpreter','none');
    end

    function IM = pickRaw()
        if S.ab == 'a', IM = S.imgs.IM_a; else, IM = S.imgs.IM_b; end
    end
    function sfc = pickSurf()
        if S.ab == 'a', sfc = S.imgs.surf_a; else, sfc = S.imgs.surf_b; end
    end
    function t = leftTitle()
        t = sprintf('Raw PIV %s  frame %s  --  %s', ...
            cfg.exp_name, pairStr(S.frames(S.idx)), upper(S.ab));
    end

    % ---- right panel ----
    function renderRightFull()
        cla(S.axR, 'reset');
        if isempty(S.data)
            text(S.axR, 0.5, 0.5, 'Click  Show ->  to compute', ...
                'Units','normalized','HorizontalAlignment','center','FontSize',12);
            axis(S.axR,'off');
            return
        end
        [fld, lim, lbl] = pickField();
        S.hImgR = imagesc(S.axR, S.data.x, S.data.z, fld);
        colormap(S.axR, S.cmap);
        colorbar(S.axR);
        clim(S.axR, lim);
        xlabel(S.axR,'x (m)'); ylabel(S.axR,'z (m)');
        title(S.axR, rightTitle(lbl), 'Interpreter','none');
        axis(S.axR, 'image');
        hold(S.axR, 'on');
        S.hQuiver = [];
        if S.quiver, drawQuiver(); end
        hold(S.axR, 'off');
    end

    function t = rightTitle(lbl)
        if isempty(S.data), t = ''; return; end
        if S.recompute
            src = sprintf('recompute  x=[%s]  z=[%s]  ov=%d%%', ...
                num2str(S.IntrWndw_x), num2str(S.IntrWndw_z), ...
                round(S.overlapFrac*100));
        else
            src = 'precomputed';
        end
        t = sprintf('%s   frame %s   (%s)', lbl, pairStr(S.data.fnum), src);
    end

    function drawQuiver()
        [X, Z] = meshgrid(S.data.x, S.data.z);
        st = max(1, round(size(S.data.dx_smooth,2) / 40));
        rs = 1:st:size(S.data.dx_smooth,1);
        cs = 1:st:size(S.data.dx_smooth,2);
        S.hQuiver = quiver(S.axR, X(rs,cs), Z(rs,cs), ...
            S.data.dx_smooth(rs,cs), S.data.dz_smooth(rs,cs), ...
            S.qscale, 'k', 'LineWidth', 0.5);
    end

    function [fld, lim, lbl] = pickField()
        switch S.field
            case 'dx_raw',    fld = S.data.dx_raw;    lbl = 'dx raw (px/frame)';
            case 'dx_smooth', fld = S.data.dx_smooth; lbl = 'dx smooth (px/frame)';
            case 'dz_raw',    fld = S.data.dz_raw;    lbl = 'dz raw (px/frame)';
            case 'dz_smooth', fld = S.data.dz_smooth; lbl = 'dz smooth (px/frame)';
            case 'u',    fld = S.data.u; lim = cfg.clim_u; lbl = 'u (m/s)'; return
            case 'w',    fld = S.data.w; lim = cfg.clim_w; lbl = 'w (m/s)'; return
            case 'dcor', fld = S.data.dcor; lbl = 'dcor (peak ratio)';
            otherwise,   fld = S.data.dx_raw; lbl = 'dx raw (px/frame)';
        end
        vals = fld(isfinite(fld) & fld ~= 0);
        if numel(vals) > 10
            lim = prctile(vals, [2 98]);
        else
            lim = [-1 1];
        end
        if lim(1) >= lim(2), lim = mean(lim) + [-1 1]; end
    end

    function updateFrameLabel()
        set(S.uiFrameLbl, 'String', ...
            sprintf('Frame %s  (%d / %d)', ...
                pairStr(S.frames(S.idx)), S.frames(S.idx), S.frames(end)));
        set(S.uiSlider,    'Value',  S.idx);
        set(S.uiFrameEdit, 'String', sprintf('%d', S.frames(S.idx)));
    end

    %% ---- navigation callbacks ----
    function setIdx(newIdx)
        newIdx = max(1, min(numel(S.frames), newIdx));
        if newIdx == S.idx, set(S.uiSlider,'Value',S.idx); return; end
        if ~isempty(S.roi) && isvalid(S.roi), S.roiPos = S.roi.Position; end
        S.idx = newIdx;
        guidata(S.fig, S);
        loadImgAndShowLeft();
    end
    function stepFrame(delta), setIdx(S.idx + delta); end

    function onFrameEdit(h,~)
        v = str2double(get(h,'String'));
        j = find(S.frames == round(v), 1);
        if isempty(j), updateFrameLabel(); else, setIdx(j); end
    end

    function toggleAB()
        if S.ab == 'a', S.ab = 'b'; else, S.ab = 'a'; end
        set(S.uiAbLbl, 'String', sprintf('Viewing: %s', upper(S.ab)));
        renderLeftBlink();
        guidata(S.fig, S);
    end

    function onIW(h,~)
        idx = get(h,'Value');
        vx = dispIW_x();  vz = dispIW_z();
        S.iw_x = vx(idx);
        S.iw_z = vz(idx);
        if ~isempty(S.roi) && isvalid(S.roi)
            p = S.roi.Position;
            cx = p(1) + p(3)/2;  cy = p(2) + p(4)/2;
            S.roi.Position = [cx - S.iw_x/2, cy - S.iw_z/2, S.iw_x, S.iw_z];
            S.roi.Label    = sprintf('IW %dx%d px', S.iw_x, S.iw_z);
            S.roiPos       = S.roi.Position;
        end
        guidata(S.fig, S);
    end
    function onROIMoved(src)
        S.roiPos = src.Position;
        guidata(S.fig, S);
    end

    function onDisp(~,~)
        xa = str2double(get(S.uiXa,'String'));
        ya = str2double(get(S.uiYa,'String'));
        xb = str2double(get(S.uiXb,'String'));
        yb = str2double(get(S.uiYb,'String'));
        if any(isnan([xa ya xb yb]))
            set(S.uiDispLbl,'String','--');
            S.ptA = [NaN NaN]; S.ptB = [NaN NaN];
        else
            ddx = xb - xa; ddz = yb - ya;
            mag = sqrt(ddx^2 + ddz^2);
            set(S.uiDispLbl,'String', ...
                sprintf('%.2f px  (dx=%.1f  dz=%.1f)', mag, ddx, ddz));
            S.ptA = [xa ya]; S.ptB = [xb yb];
        end
        if ~isempty(S.hPtA) && isvalid(S.hPtA)
            set(S.hPtA, 'XData', S.ptA(1), 'YData', S.ptA(2));
            set(S.hPtB, 'XData', S.ptB(1), 'YData', S.ptB(2));
        end
        guidata(S.fig, S);
    end

    function onField(h,~)
        S.field = cfg.fieldKeys{get(h,'Value')};
        if ~isempty(S.data), renderRightFull(); end
        guidata(S.fig, S);
    end

    function onQuiver(h,~)
        S.quiver = logical(get(h,'Value'));
        if isempty(S.data), guidata(S.fig,S); return; end
        if S.quiver
            hold(S.axR,'on'); drawQuiver(); hold(S.axR,'off');
        elseif ~isempty(S.hQuiver) && isvalid(S.hQuiver)
            delete(S.hQuiver); S.hQuiver = [];
        end
        guidata(S.fig, S);
    end

    function onQscale(h,~)
        S.qscale = get(h,'Value');
        set(S.uiQscaleLbl,'String',sprintf('%.1f', S.qscale));
        if S.quiver && ~isempty(S.data)
            if ~isempty(S.hQuiver) && isvalid(S.hQuiver), delete(S.hQuiver); end
            hold(S.axR,'on'); drawQuiver(); hold(S.axR,'off');
        end
        guidata(S.fig, S);
    end

    %% ---- Show button: load/compute velocity + update right panel ----
    function onShow(~,~)
        fnum = S.frames(S.idx);
        set(S.fig,'Pointer','watch');
        if ~isKey(S.rawCache, cacheKey(fnum))
            if S.recompute
                set(S.uiFrameLbl,'String', ...
                    sprintf('Recomputing frame %s ...', pairStr(fnum)));
            else
                set(S.uiFrameLbl,'String', ...
                    sprintf('Loading velocity %s ...', pairStr(fnum)));
            end
            drawnow;
        end
        raw    = loadRaw(fnum);
        S.data = deriveFields(raw);
        set(S.fig,'Pointer','arrow');
        updateFrameLabel();
        renderRightFull();
        guidata(S.fig, S);
    end

    %% ---- velocity source + pyramid callbacks ----
    function onSource(h,~)
        isRect = (get(h,'Value') == 2);
        S.recompute = isRect;
        set(S.uiIntrWndw_x, 'Enable', enbStr(isRect));
        set(S.uiIntrWndw_z, 'Enable', enbStr(isRect));
        set(S.uiOverlap,    'Enable', enbStr(isRect));
        % Clear velocity cache on source switch; keep image cache.
        S.rawCache = containers.Map('KeyType','char','ValueType','any');
        refreshIWPopup();
        set(S.uiGsLbl, 'String', gsLabel());
        S.roiPos = [];
        guidata(S.fig, S);
        loadImgAndShowLeft();
    end

    function s = enbStr(b)
        if b, s = 'on'; else, s = 'off'; end
    end

    function onIntrWndw_x(h,~)
        v = str2num(get(h,'String')); %#ok<ST2NM>
        v = round(v(:)).';
        if isempty(v) || any(v < 2) || any(mod(v,2)~=0) || any(diff(v)>0)
            set(h,'String', num2str(S.IntrWndw_x));
            warndlg('IW_x: even, non-increasing (e.g. 64 32 16 8).','Invalid IW_x');
            return
        end
        if numel(v) ~= numel(S.IntrWndw_z)
            set(h,'String', num2str(S.IntrWndw_x));
            warndlg(sprintf('IW_x must have %d levels (same as IW_z).', ...
                numel(S.IntrWndw_z)),'Level mismatch');
            return
        end
        S.IntrWndw_x = v;
        S.GrdSpc_x   = max(2, 2*round(v*(1-S.overlapFrac)/2));
        set(S.uiGsLbl,'String', gsLabel());
        refreshIWPopup();
        guidata(S.fig, S);
    end

    function onIntrWndw_z(h,~)
        v = str2num(get(h,'String')); %#ok<ST2NM>
        v = round(v(:)).';
        if isempty(v) || any(v < 2) || any(mod(v,2)~=0) || any(diff(v)>0)
            set(h,'String', num2str(S.IntrWndw_z));
            warndlg('IW_z: even, non-increasing (e.g. 16 8 8 4).','Invalid IW_z');
            return
        end
        if numel(v) ~= numel(S.IntrWndw_x)
            set(h,'String', num2str(S.IntrWndw_z));
            warndlg(sprintf('IW_z must have %d levels (same as IW_x).', ...
                numel(S.IntrWndw_x)),'Level mismatch');
            return
        end
        S.IntrWndw_z = v;
        S.GrdSpc_z   = max(2, 2*round(v*(1-S.overlapFrac)/2));
        set(S.uiGsLbl,'String', gsLabel());
        refreshIWPopup();
        guidata(S.fig, S);
    end

    function onOverlap(h,~)
        overlapVals = [0.5, 0.25];      % 50%  |  75%
        S.overlapFrac = overlapVals(get(h,'Value'));
        S.GrdSpc_x = max(2, 2*round(S.IntrWndw_x*(1-S.overlapFrac)/2));
        S.GrdSpc_z = max(2, 2*round(S.IntrWndw_z*(1-S.overlapFrac)/2));
        set(S.uiGsLbl,'String', gsLabel());
        guidata(S.fig, S);
    end

    function refreshIWPopup()
        str = iwStrList();
        vx  = dispIW_x();
        vz  = dispIW_z();
        nL  = numel(str);
        set(S.uiIW, 'String', str, 'Value', nL);
        S.iw_x = vx(nL);
        S.iw_z = vz(nL);
    end

    function onKey(~, ev)
        switch ev.Key
            case 'rightarrow', stepFrame(+1);
            case 'leftarrow',  stepFrame(-1);
            case 'space',      toggleAB();
            case 'u'
                S.field='u'; set(S.uiField,'Value',5);
                if ~isempty(S.data), renderRightFull(); end
                guidata(S.fig,S);
            case 'w'
                S.field='w'; set(S.uiField,'Value',6);
                if ~isempty(S.data), renderRightFull(); end
                guidata(S.fig,S);
        end
    end
end
