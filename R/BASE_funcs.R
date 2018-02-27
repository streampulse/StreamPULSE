#' Internal functions
#'
#' Not intended to be called directly by the user.
#'
#' Not intended to be called directly by the user.
#'
#' @keywords internal
#'
fit_BASE = function(directory, interval=900, n.iter=10000, n.burnin=5000){
    # CODE ADAPTED FROM https://github.com/dgiling/BASE/tree/master/BASE
    # set up output table dataframes
    filenames<-list.files(file.path(directory,"input"))
    seconds<-86400
    output.table<-NULL
    # output.table<-data.frame(File=filenames, GPP.mean=NA, GPP.sd=NA, ER.mean=NA, ER.sd=NA, K.mean=NA, K.sd=NA, theta.mean=NA, theta.sd=NA, A.mean=NA, A.sd=NA, p.mean=NA, p.sd=NA,
    #                          R2=NA, PPP=NA, rmse=NA, rmse.relative=NA, mrl.fraction=NA, ER.K.cor=NA, convergence.check=NA, A.Rhat=NA,
    #                          K.Rhat=NA, theta.Rhat=NA, p.Rhat=NA, R.Rhat=NA, GPP.Rhat=NA, DIC=NA, pD=NA)

    output.cols <- c("File","GPP.mean","GPP.sd","ER.mean","ER.sd","K.mean","K.sd","theta.mean","theta.sd","A.mean","A.sd","p.mean","p.sd",
        "R2","PPP","rmse","rmse.relative","mrl.fraction","ER.K.cor","convergence.check","A.Rhat","K.Rhat","theta.Rhat","p.Rhat","R.Rhat","GPP.Rhat","DIC","pD")
    write.table(as.matrix(t(output.cols)), file=file.path(directory,"output/BASE_results.csv"),row.names=F,col.names=F,sep=",") # output file name

    instant.rates<-data.frame(File=rep(filenames, each=seconds/interval), interval=1:(seconds/interval), tempC=NA, PAR=NA, K.instant=NA, GPP.instant=NA, ER.instant=NA)

    pb <- txtProgressBar(0,length(filenames),char="~",style=3) # create progress bar

    # Start file iteration
    for (fname in filenames){

        # fname <- filenames[1]
        data<-read.csv(file.path(directory,"input",fname), head=T) # read next file

        # define data vectors
        num.measurements<-nrow(data)
        tempC<-data$tempC
        DO.meas<-data$DO.meas
        PAR<-data$I
        salinity<-data$salinity
        atmo.pressure<-data$atmo.pressure

        inits <- function() list(sd=0.1)

        # different random seeds
        kern <- as.integer(runif(1000,min=1,max=10000))
        iters <- sample(kern,1)
        n.chains <- 3
        n.thin <- 10
        data.list <- list("num.measurements","interval","tempC","DO.meas","PAR","salinity","atmo.pressure")

        # define monitoring variables
        params <- c("A","R","K","K.day","p","theta","tau","ER","GPP", "NEP","sum.obs.resid","sum.ppa.resid","PPfit","DO.modelled", "P1", "P2")

        ## call jags

        # Set debug = T below to inspect each file for model convergence
        # (inspect the main parameters for convergence using bgr diagrams, history, density and autocorrelation)
        metab <- NULL
        metab <- do.call(jags.parallel,
            list(data=data.list, inits=inits, parameters.to.save=params, model.file = "BASE_metab_model_v2.2.txt",
                n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                n.thin = n.thin, n.cluster= n.chains, DIC = TRUE,
                working.directory = file.path(directory), jags.seed = 123, digits=5))

        # print(metab, digits=2) # to inspect results of last "metab" estimate

        ## diagnostic summaries
        # Rhat (srf) test
        srf<- metab$BUGSoutput$summary[,8]
        Rhat.test <- NULL
        Rhat.test <- ifelse(any(srf>1.1, na.rm=T)==TRUE,"Check convergence", "Fine")

        # autocorr test
        metab.mcmc<-as.mcmc(metab)
        ac.lag1 <- autocorr.diag(metab.mcmc, lags = 1)
        auto.corr.test <- NULL
        auto.corr.test <- ifelse(any(abs(ac.lag1)>0.1, na.rm=T)==TRUE,"Check ac", "ac OK")

        PPP <- metab$BUGSoutput$summary["PPfit","mean"] # posterior predictive p-value # updated

        DO.mod.means <- metab$BUGSoutput$mean$DO.modelled
        DO.mod.sd <- metab$BUGSoutput$sd$DO.modelled

        R2 <- cor(DO.mod.means,DO.meas)^2
        rmse <- sqrt(sum((metab$BUGSoutput$mean$DO.modelled-DO.meas)^2)/length(DO.meas))
        post.mean.dev <- metab$BUGSoutput$mean$deviance
        pD <- metab$BUGSoutput$pD
        DIC <- metab$BUGSoutput$DIC

        DO.lag<-DO.meas[2:length(DO.meas)]-DO.meas[1:(length(DO.meas)-1)]
        ptpvar <- sqrt((sum((DO.lag)^2)/(length(DO.meas)-1))) # point to point variation
        rmse.relative <- rmse / ptpvar

        diff<-metab$BUGSoutput$mean$DO.modelled-DO.meas
        mrl.max<-max(rle(sign(as.vector(diff)))$lengths)
        mrl.fraction<-max(rle(sign(as.vector(diff)))$lengths)/length(DO.meas) # prop of largest run

        ER.K.cor <- cor(metab$BUGSoutput$sims.list$ER,metab$BUGSoutput$sims.list$K) # plot(metab$sims.list$ER ~ metab$sims.list$K)

        # insert results to table and write table
        result <- c(fname, metab$BUGSoutput$mean$GPP, metab$BUGSoutput$sd$GPP, metab$BUGSoutput$mean$ER, metab$BUGSoutput$sd$ER, metab$BUGSoutput$mean$K.day,
            metab$BUGSoutput$sd$K.day,  metab$BUGSoutput$mean$theta, metab$BUGSoutput$sd$theta, metab$BUGSoutput$mean$A, metab$BUGSoutput$sd$A, metab$BUGSoutput$mean$p, metab$BUGSoutput$sd$p,
            R2, PPP, rmse, rmse.relative, mrl.fraction, ER.K.cor, Rhat.test, metab$BUGSoutput$summary["A",8] , metab$BUGSoutput$summary["K",8],
            metab$BUGSoutput$summary["theta",8], metab$BUGSoutput$summary["p",8], metab$BUGSoutput$summary["R",8], metab$BUGSoutput$summary["GPP",8],  DIC, pD)

        # row <- which(output.table$File==fname)
        # output.table[row,]<-result
        # write.csv(output.table, file=file.path(directory,"BASE/output/BASE_results.csv"),row.names=F) # output file name
        write.table(as.matrix(t(result)), file=file.path(directory,"output/BASE_results.csv"),row.names=F,col.names=F,append=T,sep=",")

        # insert results to instantaneous table
        rows <- which(instant.rates$File==fname)
        instant.rates$tempC[rows] <- tempC
        instant.rates$PAR[rows] <- PAR
        instant.rates$K.instant[rows] <- metab$BUGSoutput$mean$K * 1.0241^(tempC-mean(tempC))
        instant.rates$ER.instant[rows] <- metab$BUGSoutput$mean$R * metab$BUGSoutput$mean$theta^(tempC-mean(tempC))
        instant.rates$GPP.instant[rows] <- metab$BUGSoutput$mean$A * PAR^(metab$BUGSoutput$mean$p)
        write.csv(instant.rates, file=file.path(directory,"output/instantaneous rates/instantaneous_rates.csv")) # output file name

        # diagnostic multi-plot
        jpeg(file=file.path(directory,"output/validation plots", paste0(fname, ".jpg")), width=1200, height=1200, pointsize=30)
        # fit plot
        R2jags::traceplot(metab, varname=c('A','p','R','K.day','theta'), ask=FALSE, mfrow=c(3,3), mar=c(2,2,0,8), new=FALSE)
        plot(1:num.measurements,DO.mod.means, type="l",lwd=2, ylim=c(min(DO.mod.means-DO.mod.sd)-0.5,max(DO.mod.means+DO.mod.sd)+0.5), xlab="Timestep")
        points(1:num.measurements,DO.meas,pch=1,xlab="Timestep", col="grey60", cex=0.75)
        points(1:num.measurements,DO.mod.means+DO.mod.sd, type="l", lty=2)
        points(1:num.measurements,DO.mod.means-DO.mod.sd, type="l", lty=2)
        plot(1:num.measurements,tempC,pch=1,xlab="Timestep" , typ='p')
        plot(1:num.measurements,PAR,pch=1,xlab="Timestep" , typ='p')
        graphics.off()

        setTxtProgressBar(pb, which(filenames==fname)) # update progress bar

    }
    close(pb) # close progress bar

}
