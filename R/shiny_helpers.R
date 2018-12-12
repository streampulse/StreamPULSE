#' Internal functions
#'
#' Not intended to be called directly by the user.
#'
#' Not intended to be called directly by the user.
#'
#' @keywords internal
#'
processing_func = function (ts, st, en) {
    gpp = ts$GPP; gppup = ts$GPP.upper; gpplo = ts$GPP.lower
    er = ts$ER; erup = ts$ER.upper; erlo = ts$ER.lower

    ts_full = as.data.frame(ts)

    ts_full$Year = as.numeric(format(ts_full$date, "%Y"))
    ts_full$DOY = as.numeric(format(ts_full$date, "%j"))
    ts_full$NPP = ts_full$GPP + ts_full$ER

    ts_full = ts_full[-c(1,nrow(ts_full)), -c(1,8,9,10)]
    ts_full = ts_full[ts_full$DOY > st & ts_full$DOY < en,]

    return(ts_full)
}

season_ts_func = function (ts_full, fit_daily, st, en, overlay=NULL){

    ts_full = ts_full[, colnames(ts_full) != 'Year']

    avg_trajectory = aggregate(ts_full, by=list(ts_full$DOY),
        FUN=mean, na.rm=TRUE)

    gpp = avg_trajectory$GPP
    gppup = avg_trajectory$GPP.upper; gpplo = avg_trajectory$GPP.lower
    er = avg_trajectory$ER
    erup = avg_trajectory$ER.upper; erlo = avg_trajectory$ER.lower
    doy = avg_trajectory$DOY

    sd_trajectory = aggregate(ts_full, by=list(ts_full$DOY),
        FUN=sd, na.rm=TRUE)

    #get bounds
    llim = min(c(gpplo, erlo), na.rm=TRUE)
    ulim = max(c(gppup, erup), na.rm=TRUE)
    maxmin_day = range(doy, na.rm=TRUE)

    plot(doy, avg_trajectory$GPP, type="l", col="red", xlab='', las=0,
        ylab='', xaxs='i', yaxs='i',
        ylim=c(llim, ulim), lwd=2, xaxt='n', bty='l',
        xlim=c(max(st, maxmin_day[1]), min(en, maxmin_day[2])))
    mtext(expression(paste("O"[2] * " (gm"^"-2" * " d"^"-1" * ')')), side=2,
        line=2.5, font=2)

    #split time and DO series into NA-less chunks for plotting polygons
    ff = data.frame(doy=doy, gpplo=gpplo, gppup=gppup, erlo=erlo, erup=erup)
    if(! is.null(overlay) && overlay == 'mean daily K600'){
        ff$Kup = fit_daily$K600_daily_97.5pct
        ff$Klo = fit_daily$K600_daily_2.5pct
    }
    rl = rle(is.na(ff$gpplo))
    vv = !rl$values
    chunkfac = rep(cumsum(vv), rl$lengths)
    chunkfac[chunkfac == 0] = 1
    chunks = split(ff, chunkfac)
    noNAchunks = lapply(chunks, function(x) x[!is.na(x$gpplo),] )

    for(i in 1:length(noNAchunks)){
        polygon(x=c(noNAchunks[[i]]$doy, rev(noNAchunks[[i]]$doy)),
            y=c(noNAchunks[[i]]$gpplo, rev(noNAchunks[[i]]$gppup)),
            col=adjustcolor('red', alpha.f=0.3), border=NA)
        polygon(x=c(noNAchunks[[i]]$doy, rev(noNAchunks[[i]]$doy)),
            y=c(noNAchunks[[i]]$erlo, rev(noNAchunks[[i]]$erup)),
            col=adjustcolor('blue', alpha.f=0.3), border=NA)
    }

    lines(doy, avg_trajectory$ER, col="blue", lwd=2)
    abline(h=0, lty=3, col='gray50')

    #overlay user selected variable
    if(! is.null(overlay) && overlay == 'mean daily K600'){
        par(new=TRUE)
        llim2 = min(ff$Klo, na.rm=TRUE)
        ulim2 = max(ff$Kup, na.rm=TRUE)
        plot(doy, fit_daily$K600_daily_mean,
            col='orange',
            type='l', xlab='', las=0, ylab='', xaxs='i', yaxs='i',
            lwd=2, xaxt='n', bty='u', yaxt='n', ylim=c(llim2, ulim2),
            xlim=c(max(st, maxmin_day[1]), min(en, maxmin_day[2])))
        axis(4)#, col.axis='orange')
        mtext(expression(paste('K600 (d'^'-1' * ')')), side=4,
            line=2.5, font=2)#, col='orange')
        for(i in 1:length(noNAchunks)){
            polygon(x=c(noNAchunks[[i]]$doy, rev(noNAchunks[[i]]$doy)),
                y=c(noNAchunks[[i]]$Klo, rev(noNAchunks[[i]]$Kup)),
                col=adjustcolor('orange', alpha.f=0.3), border=NA)
        }
    }
}

metab_legend = function(show_K600=FALSE){

    # par(rep(0,4), oma=rep(0,4))
    par(mar=c(0,4,0,1), oma=rep(0,4))
    plot(1,1, axes=FALSE, type='n', xlab='', ylab='', bty='o')
    if(show_K600 == FALSE){
        legend("bottomleft", ncol=2, xpd=FALSE,
            legend=c("GPP", "ER"), bty="n", lty=1,
            lwd=2, col=c("red", "blue"),
            x.intersp=c(.5,.5))
        legend('bottom', horiz=TRUE, seg.len=1,
            bty="n", lty=1, legend=c('95% CIs',''),
            col=c(adjustcolor('red', alpha.f=0.3),
                adjustcolor('blue', alpha.f=0.3)),
            lwd=6)
    } else {
        legend("bottomleft", ncol=3, xpd=FALSE,
            legend=c("GPP", "ER", 'K600'), bty="n", lty=1,
            lwd=2, col=c("red", "blue", 'orange'),
            x.intersp=c(.5,.5,.5))
        legend("bottom", ncol=3, xpd=FALSE,
            legend=c('', '', '95% CIs'), bty='n', lty=1, lwd=6,
            col=c(adjustcolor('red', alpha.f=0.3),
                adjustcolor('blue', alpha.f=0.3),
                adjustcolor('orange', alpha.f=0.3)),
            x.intersp=c(0,0,0), text.width=c(0,0,0))
    }
}

kernel_func = function(ts_full, main){

    kernel = ks::kde(na.omit(ts_full[, c('GPP', 'ER')]))
    k_lim = max(abs(c(min(ts_full$ER, na.rm=TRUE),
        max(ts_full$GPP, na.rm=TRUE))))
    plot(kernel, xlab='', las=1, xaxt='n', ylab='', yaxt='n',
        ylim=c(-k_lim, 0), xlim=c(0, k_lim), display='filled.contour',
        col=c(NA, "purple1", "purple3", "purple4"))
    axis(1, tcl=-0.2, padj=-1)
    axis(2, tcl=-0.2, hadj=0.5, las=1)
    mtext(expression(paste("GPP (gm"^"-2" * " d"^"-1" * ")")),
        1, line=1.8)
    mtext(expression(paste("ER (gm"^"-2" * " d"^"-1" * ")")),
        2, line=2)
    abline(0, -1, col='black', lty=3)
}

kernel_legend = function(){
    par(mar = rep(0,4), oma = rep(0,4))
    plot(1,1, axes=FALSE, type='n', xlab='', ylab='', bty='o')
    legend("bottomright", c("75%", "50%", "25%"), bty="o", bg='white',
        lty=c(1,1,1), lwd=4, col=c("purple1", "purple3", "purple4"),
        seg.len=1, box.col='transparent', horiz=TRUE)
}

O2_plot = function(mod_out, st, en, brush, click, overlay='None',
    xformat='DOY', varmap){
    # st=0; en=366
    # brush = list(xmin=1460617508, xmax=1464058124, ymin=8.816155, ymax=14.45195)

    #convert POSIX time to DOY and UNIX time
    DOY = as.numeric(gsub('^0+', '', strftime(mod_out$data$solar.time,
        format="%j")))
    ustamp = as.numeric(as.POSIXct(mod_out$data$solar.time, tz='UTC'))

    #replace initial DOYs of 365 or 366 (solar date in previous calendar year) with 1
    if(DOY[1] %in% 365:366){
        DOY[DOY %in% 365:366 & 1:length(DOY) < length(DOY)/2] = 1
    }

    #get bounds
    xmin_ind = match(st, DOY)
    if(is.na(xmin_ind)) xmin_ind = 1
    xmin = ustamp[xmin_ind]

    xmax_ind = length(DOY) - match(en, rev(DOY)) + 1
    if(is.na(xmax_ind)) xmax_ind = nrow(mod_out$data)
    xmax = ustamp[xmax_ind]

    #overlay additional series if selected
    if(overlay != 'None'){
        cleaned_varnames = sapply(varmap, function(x) x[[1]])
        overlay = names(varmap)[which(cleaned_varnames == overlay)]
        slice = mod_out$data[xmin_ind:xmax_ind, c('DO.obs', 'DO.mod', overlay)]
    } else {
        slice = mod_out$data[xmin_ind:xmax_ind, c('DO.obs', 'DO.mod')]
    }
    yrng = range(c(slice$DO.obs, slice$DO.mod), na.rm=TRUE)

    #window, series, axis labels
    plot(ustamp, mod_out$data$DO.obs, xaxt='n', las=0,
        type='n', xlab='', ylab='', bty='l', ylim=c(yrng[1], yrng[2]),
        xaxs='i', yaxs='i', xlim=c(xmin, xmax))

    #split time and DO series into NA-less chunks for plotting polygons
    ff = data.frame(ustamp=ustamp, DO=mod_out$data$DO.mod,
        zero=rep(0, length(mod_out$data$DO.mod)))
    rl = rle(is.na(ff$DO))
    vv = !rl$values
    chunkfac = rep(cumsum(vv), rl$lengths)
    chunkfac[chunkfac == 0] = 1
    chunks = split(ff, chunkfac)
    noNAchunks = lapply(chunks, function(x) x[!is.na(x$DO),] )

    for(i in 1:length(noNAchunks)){
        polygon(x=c(noNAchunks[[i]]$ustamp, rev(noNAchunks[[i]]$ustamp)),
            y=c(noNAchunks[[i]]$DO, rep(0, nrow(noNAchunks[[i]]))),
            col='gray75', border='gray75')
    }

    mtext(expression(paste('DO (mgL'^'-1' * ')')), 2, font=1, line=2.5)
    lines(ustamp, mod_out$data$DO.obs, col='cyan4')

    #get seq of 10 UNIX timestamps and use corresponding DOYs/dates as ticks
    tcs = seq(xmin, xmax, length.out=10)
    near_ind = findInterval(tcs, ustamp)

    if(xformat == 'DOY'){
        axis(1, ustamp[near_ind], DOY[near_ind], tcl=-0.2, padj=-1)
        mtext('DOY', 1, font=1, line=1.5)
    } else { #xformat == 'Date'
        date = as.Date(gsub('^0+', '', strftime(mod_out$data$solar.time,
            format="%Y-%m-%d")))
        if(DOY[1] %in% 365:366){
            date[DOY %in% 365:366 & 1:length(DOY) < length(DOY)/2] = NA
        }
        axis(1, ustamp[near_ind], date[near_ind], tcl=-0.2, padj=-1)
        mtext('Date', 1, font=1, line=1.5)
    }

    # overlay user selected variable
    if(overlay != 'None'){
        par(new=TRUE)
        yrng2 = range(slice[,overlay], na.rm=TRUE)
        plot(ustamp, mod_out$data[,overlay], xaxt='n', las=0, yaxt='n',
            type='l', xlab='', ylab='', bty='u', ylim=c(yrng2[1], yrng2[2]),
            xaxs='i', yaxs='i', xlim=c(xmin, xmax), col='coral4')
        axis(4)
        mtext(varmap[[overlay]][[2]], side=4, line=2.5)
    }

    #highlight brushed points
    hl_log_ind = which(mod_out$data$DO.mod < brush$ymax &
            mod_out$data$DO.mod > brush$ymin &
            ustamp < brush$xmax & ustamp > brush$xmin)
    hl_x = ustamp[hl_log_ind]
    hl_y = mod_out$data$DO.mod[hl_log_ind]
    points(hl_x, hl_y, col='goldenrod1', cex=0.3, pch=20)

}

O2_legend = function(overlay, varmap){
    par(mar=c(0,4,0,1), oma=rep(0,4))
    plot(1,1, axes=FALSE, type='n', xlab='', ylab='', bty='o')
    if(overlay == 'None'){
        legend(x='bottomleft', legend=c('Pred', 'Obs'), bg='white',
            cex=0.8, col=c('gray75', 'cyan4'), lty=1, bty='o', horiz=TRUE,
            lwd=c(6,2), box.col='transparent')
    } else {
        cleaned_varnames = sapply(varmap, function(x) x[[1]])
        overlay = names(varmap)[which(cleaned_varnames == overlay)]
        legend(x='bottomleft', legend=c('Pred DO', 'Obs DO',
            varmap[[overlay]][[1]]),
            bg='white', cex=0.8, col=c('gray75', 'cyan4', 'coral4'),
            lty=1, bty='o', horiz=TRUE,
            lwd=c(6,2,2), box.col='transparent')
    }
}

KvER_plot = function(mod_out, slice, click=NULL){

    mod = lm(slice$ER_mean ~ slice$K600_daily_mean)
    R2 = sprintf('%1.2f', summary(mod)$adj.r.squared)
    plot(slice$K600_daily_mean, slice$ER_mean,
        col='darkgreen', ylab='', xlab='Daily mean K600',
        bty='l', font.lab=1, cex.axis=0.8, las=1)
    mtext(expression(paste("Daily mean ER (gm"^"-2" ~ ")")), side=2, line=2.5)
    mtext(bquote('Adj.' ~ R^2 * ':' ~ .(R2)), side=3, line=0, adj=1,
        cex=0.8, col='gray50')
    abline(mod, lty=2, col='gray50', lwd=2)

    #highlight point and display date on click
    if(! is.null(click) && ! is.null(click$x)){
        xmax = max(slice$K600_daily_mean, na.rm=TRUE)
        xmin = min(slice$K600_daily_mean, na.rm=TRUE)
        ymax = max(slice$ER_mean, na.rm=TRUE)
        ymin = min(slice$ER_mean, na.rm=TRUE)
        xrng = xmax - xmin
        yrng = ymax - ymin

        click_ind = which(slice$ER_mean < click$y + 0.01 * yrng &
                slice$ER_mean > click$y - 0.01 * yrng &
                slice$K600_daily_mean < click$x + 0.01 * xrng &
                slice$K600_daily_mean > click$x - 0.01 * xrng)[1]
        click_x = slice$K600_daily_mean[click_ind]
        click_y = slice$ER_mean[click_ind]
        points(click_x, click_y, col='goldenrod1', pch=19, cex=2)
        x_prop = scales::rescale(click_x, c(0, 1), c(xmin, xmax))
        if(! is.na(x_prop) && x_prop > 0.5){
            text(click_x, click_y, '-', pos=2, font=2, col='white', cex=9)
            text(click_x, click_y, '- ', pos=2, font=2, col='white', cex=9)
            text(click_x, click_y, paste0(slice$date[click_ind], ' '), pos=2, font=2)
        } else {
            text(click_x, click_y, '-', pos=4, font=2, col='white', cex=9)
            text(click_x, click_y, ' -', pos=4, font=2, col='white', cex=9)
            text(click_x, click_y, '--', pos=4, font=2, col='white', cex=9)
            text(click_x, click_y, paste0(' ', slice$date[click_ind]), pos=4, font=2)
        }
    }
}

KvGPP_plot = function(mod_out, slice, click=NULL){

    plot(slice$K600_daily_mean, slice$GPP_mean,
        col='darkblue', ylab='', xlab='Daily mean K600',
        bty='l', font.lab=1, cex.axis=0.8, las=1)
    mtext(expression(paste("Daily mean GPP (gm"^"-2" ~ ")")), side=2, line=2.5)

    #highlight point and display date on click
    if(! is.null(click) && ! is.null(click$x)){
        xmax = max(slice$K600_daily_mean, na.rm=TRUE)
        xmin = min(slice$K600_daily_mean, na.rm=TRUE)
        ymax = max(slice$GPP_mean, na.rm=TRUE)
        ymin = min(slice$GPP_mean, na.rm=TRUE)
        xrng = xmax - xmin
        yrng = ymax - ymin

        click_ind = which(slice$GPP_mean < click$y + 0.01 * yrng &
                slice$GPP_mean > click$y - 0.01 * yrng &
                slice$K600_daily_mean < click$x + 0.01 * xrng &
                slice$K600_daily_mean > click$x - 0.01 * xrng)[1]
        click_x = slice$K600_daily_mean[click_ind]
        click_y = slice$GPP_mean[click_ind]
        points(click_x, click_y, col='goldenrod1', pch=19, cex=2)
        x_prop = scales::rescale(click_x, c(0, 1), c(xmin, xmax))
        if(! is.na(x_prop) && x_prop > 0.5){
            text(click_x, click_y, '-', pos=2, font=2, col='white', cex=9)
            text(click_x, click_y, '- ', pos=2, font=2, col='white', cex=9)
            text(click_x, click_y, paste0(slice$date[click_ind], ' '), pos=2, font=2)
        } else {
            text(click_x, click_y, '-', pos=4, font=2, col='white', cex=9)
            text(click_x, click_y, ' -', pos=4, font=2, col='white', cex=9)
            text(click_x, click_y, paste0(' ', slice$date[click_ind]), pos=4, font=2)
        }
    }
}

KvQ_plot = function(mod_out, slicex, slicey, click=NULL){

    nodes = mod_out$fit$KQ_binned$lnK600_lnQ_nodes_mean
    log_Q = log(slicex$discharge.daily)
    xminplot = min(c(log_Q, nodes), na.rm=TRUE)
    xmaxplot = max(c(log_Q, nodes), na.rm=TRUE)
    KQmod = lm(slicey$K600_daily_mean ~ log_Q)
    R2 = sprintf('%1.2f', summary(KQmod)$adj.r.squared)
    plot(log_Q, slicey$K600_daily_mean, xlim=c(xminplot, xmaxplot),
        col='purple4', xlab='Log daily mean Q (cms)', ylab='Daily mean K600',
        bty='l', font.lab=1, cex.axis=0.8, las=1)
    abline(v=nodes, lty=2, col='darkred')
    mtext('log Q node centers', side=3, line=0, adj=1, cex=0.8, col='darkred')
    mtext(bquote('Adj.' ~ R^2 * ':' ~ .(R2)), side=3, line=0, adj=0,
        cex=0.8, col='gray50')
    abline(KQmod, lty=2, col='gray50', lwd=2)
    # legend('topright', legend='log Q node centers', bty='n', lty=2, col='gray')

    #highlight point and display date on click
    if(! is.null(click) && ! is.null(click$x)){
        xmax = max(log_Q, na.rm=TRUE)
        xmin = min(log_Q, na.rm=TRUE)
        ymax = max(slicey$K600_daily_mean, na.rm=TRUE)
        ymin = min(slicey$K600_daily_mean, na.rm=TRUE)
        xrng = xmax - xmin
        yrng = ymax - ymin

        click_ind_x = which(log_Q < click$x + 0.01 * xrng & log_Q >
                click$x - 0.01 * xrng)[1]
        click_ind_y = which(slicey$K600_daily_mean < click$y + 0.01 * yrng &
                slicey$K600_daily_mean > click$y - 0.01 * yrng)[1]
        click_x = log_Q[click_ind_x]
        click_y = slicey$K600_daily_mean[click_ind_y]
        points(click_x, click_y, col='goldenrod1', pch=19, cex=2)
        x_prop = scales::rescale(click_x, c(0, 1), c(xmin, xmax))
        if(! is.na(x_prop) && x_prop > 0.5){
            text(click_x, click_y, '-', pos=2, font=2, col='white', cex=9)
            text(click_x, click_y, '- ', pos=2, font=2, col='white', cex=9)
            text(click_x, click_y, paste0(slicey$date[click_ind_y], ' '), pos=2, font=2)
        } else {
            text(click_x, click_y, '-', pos=4, font=2, col='white', cex=9)
            text(click_x, click_y, ' -', pos=4, font=2, col='white', cex=9)
            text(click_x, click_y, paste0(' ', slicey$date[click_ind_y]), pos=4, font=2)
        }
    }
}

QvKres_plot = function(mod_out, slicex, slicey, click=NULL){

    log_Q = log(slicex$discharge.daily)

    #get K residuals (based on K-Q linear relationship)
    KQmod = lm(slicey$K600_daily_mean ~ log_Q, na.action=na.exclude)
    KQrelat = fitted(KQmod)
    # if(length(KQrelat) != length(slicey$K600_daily_mean)){
    #     if(is.na(slicey$K600_daily_mean[1])){
    #         KQrelat = c(NA, KQrelat)
    #     } else {
    #         KQrelat = c(KQrelat, NA)
    #     }
    # }
    # if(length(KQrelat) != length(slicey$K600_daily_mean)){
    #     KQrelat = c(KQrelat, NA)
    # }

    resid = slicey$K600_daily_mean - KQrelat

    plot(log_Q, resid,
        col='purple4', xlab='Log daily mean Q (cms)',
        ylab='Daily mean K600 residuals*',
        bty='l', font.lab=1, cex.axis=0.8, las=1)

    #highlight point and display date on click
    if(! is.null(click) && ! is.null(click$x)){
        xmax = max(log_Q, na.rm=TRUE)
        xmin = min(log_Q, na.rm=TRUE)
        ymax = max(resid, na.rm=TRUE)
        ymin = min(resid, na.rm=TRUE)
        xrng = xmax - xmin
        yrng = ymax - ymin

        click_ind_x = which(log_Q < click$x + 0.01 * xrng & log_Q >
                click$x - 0.01 * xrng)[1]
        click_ind_y = which(resid < click$y + 0.01 * yrng &
                resid > click$y - 0.01 * yrng)[1]
        click_x = log_Q[click_ind_x]
        click_y = resid[click_ind_y]
        points(click_x, click_y, col='goldenrod1', pch=19, cex=2)
        x_prop = scales::rescale(click_x, c(0, 1), c(xmin, xmax))
        if(! is.na(x_prop) && x_prop > 0.5){
            text(click_x, click_y, '-', pos=2, font=2, col='white', cex=9)
            text(click_x, click_y, '- ', pos=2, font=2, col='white', cex=9)
            text(click_x, click_y, paste0(slicey$date[click_ind_y], ' '), pos=2, font=2)
        } else {
            text(click_x, click_y, '-', pos=4, font=2, col='white', cex=9)
            text(click_x, click_y, ' -', pos=4, font=2, col='white', cex=9)
            text(click_x, click_y, paste0(' ', slicey$date[click_ind_y]), pos=4, font=2)
        }
    }
}
