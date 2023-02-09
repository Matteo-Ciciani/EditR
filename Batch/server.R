shinyServer(
  
  function(input, output, session) {
    ############################################################################
    # READ INPUT
    ############################################################################
    # reading in file
    output$fileInputs=renderUI({
      html_ui = " "
      for (i in 1:input$nfiles){
          html_ui <- paste0(html_ui, fileInput(paste0("file",i), label=paste0("Upload .ab1 File ", i)))
          html_ui <- paste0(html_ui, textInput(paste0("guide",i), label = paste0("Enter guide sequence ", i, " in 5'-3' orientation")))
          html_ui <- paste0(html_ui, checkboxInput(paste0("guide.rev",i), "Guide sequence is reverse complement", FALSE))
          html_ui <- paste0(html_ui, '<hr style="height:1px;border:none;color:#bbb;background-color:#bbb;" />')
      }
      HTML(html_ui)
    })
    
    input.seqReactive <- reactive({
        if(!is.null(input$file1) & !is.null(input$scramble)) {
        input_files <- lapply(c(1:input$nfiles), function(i) {eval(parse(text=paste0('input$file', i)))})
        input_seqs <- lapply(input_files, function(x) readsangerseq(x$datapath))
        #input_seqs <- readsangerseq(input_files[[1]]$datapath)
        return(c(input_seqs, readsangerseq(input$scramble$datapath)))
      } else return(validate(
        need(input$file1, "Please upload your sanger sequence files"),
        need(input$scramble, "Please upload your scramble sanger sequence file")
      ))
    })
    
    # making basecalls
    
    input.basecallsReactive <- reactive({lapply(input.seqReactive(), makeBaseCalls)})
    
    # getting the guide sequence
    output$guideSeqs <- reactive(lapply(c(1:input$nfiles), function(i) {eval(parse(text=paste0('input$guide', i)))}))
    
    guideReactive <- reactive({
        validate(
          need(input$guide1 != "", "Please enter a guide RNA sequence")
        )
        input_guides <- lapply(c(1:input$nfiles), function(i) {eval(parse(text=paste0('input$guide', i)))})
        input_guides_rev <- lapply(c(1:input$nfiles), function(i) {eval(parse(text=paste0('input$guide.rev', i)))})
        guides <- lapply(input_guides, DNAString)
        guides <- mapply(function(x,y)  if (y) reverseComplement(x) else x , guides, input_guides_rev)
        return(guides)
    })
    
    ## pvalue cutoff
    
    p.val.Reactive <- reactive({
      validate(
        need(input$pvalcutoff != "", "Please enter a p-value cutoff")
      )
      input$pvalcutoff
    })
    
    filename_out <- reactive(lapply(c(1:input$nfiles), function(i) strsplit(eval(parse(text=paste0('input$file', i, '$name'))), '\\.')[[1]][1]))
    
    ############################################################################
    # PREPROCESS SANGER INPUT
    ############################################################################
    sangsReactive <- reactive({
      input.seq <- input.seqReactive()
      input.basecalls <- input.basecallsReactive()
      #input.peakampmatrix <- lapply(input.basecalls, peakAmpMatrix)
      
      ### creating a sanger sequencing data.frame
      #sangs <- mapply(CreateSangs, input.peakampmatrix, input.basecalls)
      sangs <- lapply(input.basecalls, function (x) CreateSangs(peakAmpMatrix(x), x))
      return(sangs)
    })
    
    sangs.filtReactive <- reactive({
      sangs <- sangsReactive()
      
      if(is.na(input$trim3) & is.na(input$trim5)){
        # Remove 5' poor sequencing due to poor primer binding
        sangs.filt <- lapply(sangs, function(x) x %>% filter(index > 20))
        # removing crappy end
        peakTotAreaCutoff <- lapply(sangs.filt, function(x) mean(x$Tot.area)/10)
        #sangs.filt <- lapply(sangs.filt, function(x) x %>% filter(Tot.area > peakTotAreaCutoff))
        sangs.filt <- mapply(function(x, y) x %>% filter(Tot.area > y), sangs.filt, peakTotAreaCutoff, SIMPLIFY = FALSE)
        return(sangs.filt)
      }else{
        validate(need(!is.na(input$trim3), label = "A value to trim the 3' end "))
        validate(need(!is.na(input$trim5), label = "A value to trim the 5' end "))

        sangs.filt <- lapply(sangs, function(x) x[input$trim5:input$trim3, ])
        return(sangs.filt)
      }
    })
    
    output$trimmedrange <- renderText({
      
      if(!is.na(input$trim5) & !is.na(input$trim3)){
        paste("Trimmed input, starting at ", input$trim5, " and ending at ", input$trim3)
      }

    })

    # finding the guide coordinates
    guide.coordReactive <- reactive({
      # the guide coordinates are relative to the index of the sequence that came out of the 
      # makeBaseCalls function
      # the guide is matched to the sequence of the filtered sanger sequencing, so that it doesn't
      # match to any of the low quality regions. 
      
      sangs.filt <- sangs.filtReactive()[1:input$nfiles]
      
      # getting the sequence to match to
      filt.sequence <- lapply(sangs.filt, function(x) x$base.call %>% paste(collapse = "") %>% DNAString())
      
      # this function finds where the guide matches
      guides <- guideReactive()
      guide.match <- mapply(GetGuideMatch, guides, filt.sequence, SIMPLIFY = FALSE)
      
      # Finding the index values
      guide.coord <- mapply(function(x,y) list(start = x[y$start, "index"], end = x[y$end, "index"]), sangs.filt, guide.match, SIMPLIFY = FALSE)
      return(guide.coord)
    })
    guide.coordReactiveScramble <- reactive({
        # the guide coordinates are relative to the index of the sequence that came out of the 
        # makeBaseCalls function
        # the guide is matched to the sequence of the filtered sanger sequencing, so that it doesn't
        # match to any of the low quality regions. 
        
        sangs.filt <- sangs.filtReactive()[[input$nfiles+1]]
        
        # getting the sequence to match to
        filt.sequence <- sangs.filt$base.call %>% paste(collapse = "") %>% DNAString()
        
        # this function finds where the guide matches
        guides <- guideReactive()
        guide.match <- lapply(guides, function(x) GetGuideMatch(x, filt.sequence))
        
        # Finding the index values
        guide.coord <- lapply(guide.match, function(x) list(start = sangs.filt[x$start, "index"], end = sangs.filt[x$end, "index"]))
        return(guide.coord)
    })

    ############################################################################
    # GET MODEL NULL DIST
    ############################################################################
    
    nullparams.Reactive <- reactive({
      sangs.filt <- sangs.filtReactive()[1:input$nfiles]
      guide.coord <- guide.coordReactive()
      # getting params for the different null models
      null.m.params <- mapply(GetNullDistModel, sangs.filt, guide.coord, SIMPLIFY = FALSE)
    })
    
    nullparams.ReactiveScramble <- reactive({
        sangs.filt <- rep(sangs.filtReactive()[input$nfiles+1], input$nfiles)
        guide.coord <- guide.coordReactiveScramble()
        # getting params for the different null models
        null.m.params <- mapply(GetNullDistModel, sangs.filt, guide.coord, SIMPLIFY = FALSE)
    })
    
    # Editing reactive
    # this produces a dataframe of the guide region with the probabilities that the nonprimary base shows evidence of editing
    
    editing.Reactive <- reactive({
      # grabbing the reactive stuff
      sangs.filt <- sangs.filtReactive()[1:input$nfiles]
      sangs <- sangsReactive()[1:input$nfiles]
      guide.coord <- guide.coordReactive()
      guide <- guideReactive()
      null.m.params <- nullparams.Reactive()
      
      editing.df <- mapply(CreateEditingDF, guide.coord, guide, sangs, null.m.params, SIMPLIFY = FALSE)
      return(editing.df)
    })
    editing.ReactiveScramble <- reactive({
        # grabbing the reactive stuff
        sangs.filt <- rep(sangs.filtReactive()[input$nfiles+1], input$nfiles)
        sangs <- rep(sangsReactive()[input$nfiles+1], input$nfiles)
        guide.coord <- guide.coordReactiveScramble()
        guide <- guideReactive()
        null.m.params <- nullparams.ReactiveScramble()
        
        editing.df <- mapply(CreateEditingDF, guide.coord, guide, sangs, null.m.params, SIMPLIFY = FALSE)
        return(editing.df)
    })
    
    base.infoReactive <- list()
    base.infoReactiveScramble <- list()
    
    editPlot <- list()
    editPlotScramble <- list()
    
    base.editingtabledata <- list()
    base.editingtabledataScramble <- list()
    
    base.preocessdData <- list()
    base.preocessdDataEnv <- list()
    env <- environment()

    ############################################################################
    # MAKE OUTPUT
    ############################################################################
    
    observe({
        req(input$nfiles)
        # maybe check if files and guides are present
        lapply(1:input$nfiles, function(i) {
        
        ############################################################################
        # MAKE QC PLOTS
        ############################################################################
        output[[paste0('prefilter.totalarea',i)]] <- renderPlot({
            guide.coord <- guide.coordReactive()[[i]]
            sangs <- sangsReactive()[[i]]
            ggplot(sangs, aes(x = index, y = Tot.area)) +
                geom_rect(xmin = guide.coord$start, xmax = guide.coord$end, ymin = 0, ymax = Inf, fill = "lightgrey") +
                geom_line() +
                labs(x = "Base position",
                    y = "Total peak area at base position (A + T + C + G)",
                    title = "Unfiltered data, total peak area")
        })
        output[[paste0('prefilter.totalareaScramble',i)]] <- renderPlot({
            guide.coord <- guide.coordReactiveScramble()[[i]]
            sangs <- sangsReactive()[[input$nfiles+1]]
            ggplot(sangs, aes(x = index, y = Tot.area)) +
                geom_rect(xmin = guide.coord$start, xmax = guide.coord$end, ymin = 0, ymax = Inf, fill = "lightgrey") +
                geom_line() +
                labs(x = "Base position",
                    y = "Total peak area at base position (A + T + C + G)",
                    title = "Unfiltered data, total peak area")
        })
        
        output[[paste0('postfilter.signal.noise',i)]] <- renderPlotly({
            guide.coord <- guide.coordReactive()[[i]]
            sangs.filt <- sangs.filtReactive()[[i]]
            sangs <- sangsReactive()[[i]]
            
            ## doing some data rearrangement for plotting
            
            catch <- sangs.filt %>% dplyr::select(A.area:T.area, A.perc:T.perc, base.call, index)
            catch %<>% gather( key = base, value = value,
                A.area:T.area, A.perc:T.perc) %>%
                separate(col = base, into = c("base", "measure"))
            
            noise <- catch %>% 
                group_by(index, measure) %>% 
                filter(base != base.call) %>%
                summarize(noise = sum(value))
            signal <- catch %>%
                group_by(index, measure) %>%
                filter(base == base.call) %>% 
                summarize(signal = sum(value))
            
            summary.df <- left_join(noise,signal) %>%
                gather(key = type, value = value, noise, signal)
            
            sangs.plot <- sangs.filt %>% dplyr::select(index, base.call) 
            sangs.plot %<>% left_join(summary.df)
            
            ## now acutally making the plot
            
            sangs.plot$type <- ordered(sangs.plot$type, levels = c("signal", "noise"))
            
            max.height <- filter(sangs.plot, measure == "area") %>% select(value) %>% max
            
            guide.region <- data.frame(index = c(guide.coord$start, guide.coord$end,
                guide.coord$end, guide.coord$start), 
                value = c(0, 0, max.height, max.height))
            
            p <- sangs.plot %>% filter(measure == "area") %>%
                ggplot(aes(x = index, y = value)) +
                geom_polygon(data = guide.region, fill = "lightgrey") +
                geom_line(aes(color = type)) +
                scale_color_manual(values = c("#5e3c99", "#e66101")) +  
                labs(title = "Filtered data signal and noise total area",
                    x = "Position",
                    y = "Peak Area") + 
                theme(legend.title = element_blank())
            ggplotly(p)
            
        })
        output[[paste0('postfilter.signal.noiseScramble',i)]] <- renderPlotly({
            guide.coord <- guide.coordReactiveScramble()[[i]]
            sangs.filt <- sangs.filtReactive()[[input$nfiles+1]]
            sangs <- sangsReactive()[[input$nfiles+1]]
            
            ## doing some data rearrangement for plotting
            
            catch <- sangs.filt %>% dplyr::select(A.area:T.area, A.perc:T.perc, base.call, index)
            catch %<>% gather( key = base, value = value,
                A.area:T.area, A.perc:T.perc) %>%
                separate(col = base, into = c("base", "measure"))
            
            noise <- catch %>% 
                group_by(index, measure) %>% 
                filter(base != base.call) %>%
                summarize(noise = sum(value))
            signal <- catch %>%
                group_by(index, measure) %>%
                filter(base == base.call) %>% 
                summarize(signal = sum(value))
            
            summary.df <- left_join(noise,signal) %>%
                gather(key = type, value = value, noise, signal)
            
            sangs.plot <- sangs.filt %>% dplyr::select(index, base.call) 
            sangs.plot %<>% left_join(summary.df)
            
            ## now acutally making the plot
            
            sangs.plot$type <- ordered(sangs.plot$type, levels = c("signal", "noise"))
            
            max.height <- filter(sangs.plot, measure == "area") %>% select(value) %>% max
            
            guide.region <- data.frame(index = c(guide.coord$start, guide.coord$end,
                guide.coord$end, guide.coord$start), 
                value = c(0, 0, max.height, max.height))
            
            p <- sangs.plot %>% filter(measure == "area") %>%
                ggplot(aes(x = index, y = value)) +
                geom_polygon(data = guide.region, fill = "lightgrey") +
                geom_line(aes(color = type)) +
                scale_color_manual(values = c("#5e3c99", "#e66101")) +  
                labs(title = "Filtered data signal and noise total area",
                    x = "Position",
                    y = "Peak Area") + 
                theme(legend.title = element_blank())
            ggplotly(p)
            
        })
        
        output[[paste0('postfilter.noise.perc',i)]] <- renderPlot({
            guide.coord <- guide.coordReactive()[[i]]
            sangs.filt <- sangs.filtReactive()[[i]]
            
            ## doing some data rearrangement for plotting
            
            catch <- sangs.filt %>% dplyr::select(A.area:T.area, A.perc:T.perc, base.call, index)
            catch %<>% gather( key = base, value = value,
                A.area:T.area, A.perc:T.perc) %>%
                separate(col = base, into = c("base", "measure"))
            
            # splitting the catch dataframe into either signal and noise, and calculating the 
            # total noise area or total noise percent
            noise <- catch %>% 
                group_by(index, measure) %>% 
                filter(base != base.call) %>%
                summarize(noise = sum(value))
            signal <- catch %>%
                group_by(index, measure) %>%
                filter(base == base.call) %>% 
                summarize(signal = sum(value))
            
            signal.noise.df <- left_join(noise,signal) %>%
                gather(key = type, value = value, noise, signal)
            
            # making the plotting df
            sangs.plot <- sangs.filt %>% dplyr::select(index, base.call) 
            sangs.plot %<>% left_join(signal.noise.df)
            sangs.plot$type <- ordered(sangs.plot$type, levels = c("signal", "noise"))
            
            sangs.plot %>% filter(measure == "perc", type == "noise") %>%
                ggplot(aes(x = index, y = value)) +
                geom_area(aes(fill = type)) +
                scale_fill_manual(values = c("#e66101")) + 
                annotate("rect", xmin=guide.coord$start, xmax=guide.coord$end,
                    ymin=0, ymax=Inf, alpha=1/5, fill="black") +
                guides(fill = FALSE) + 
                labs(title = "Percent peak area noise",
                    x = "Position",
                    y = "Percent peak area noise")
        })
        output[[paste0('postfilter.noise.percScramble',i)]] <- renderPlot({
            guide.coord <- guide.coordReactiveScramble()[[i]]
            sangs.filt <- sangs.filtReactive()[[input$nfiles+1]]
            
            ## doing some data rearrangement for plotting
            
            catch <- sangs.filt %>% dplyr::select(A.area:T.area, A.perc:T.perc, base.call, index)
            catch %<>% gather( key = base, value = value,
                A.area:T.area, A.perc:T.perc) %>%
                separate(col = base, into = c("base", "measure"))
            
            # splitting the catch dataframe into either signal and noise, and calculating the 
            # total noise area or total noise percent
            noise <- catch %>% 
                group_by(index, measure) %>% 
                filter(base != base.call) %>%
                summarize(noise = sum(value))
            signal <- catch %>%
                group_by(index, measure) %>%
                filter(base == base.call) %>% 
                summarize(signal = sum(value))
            
            signal.noise.df <- left_join(noise,signal) %>%
                gather(key = type, value = value, noise, signal)
            
            # making the plotting df
            sangs.plot <- sangs.filt %>% dplyr::select(index, base.call) 
            sangs.plot %<>% left_join(signal.noise.df)
            sangs.plot$type <- ordered(sangs.plot$type, levels = c("signal", "noise"))
            
            sangs.plot %>% filter(measure == "perc", type == "noise") %>%
                ggplot(aes(x = index, y = value)) +
                geom_area(aes(fill = type)) +
                scale_fill_manual(values = c("#e66101")) + 
                annotate("rect", xmin=guide.coord$start, xmax=guide.coord$end,
                    ymin=0, ymax=Inf, alpha=1/5, fill="black") +
                guides(fill = FALSE) + 
                labs(title = "Percent peak area noise",
                    x = "Position",
                    y = "Percent peak area noise")
        })
        
        ############################################################################
        # FIT MODEL
        ############################################################################
        
        base.infoReactive[[i]] <- reactive({
            # this reactive is collecting information on each base:
            # - average percent signal
            # critical value for each base for the cutoff of significance
            
            sangs.filt <- sangs.filtReactive()[[i]]
            null.m.params <- nullparams.Reactive()[[i]]
            p.val.cutoff <- p.val.Reactive()
            
            # finding the average percent signal for each base
            avg.base <- sangs.filt %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc) %>%
                separate(col = focal.base, into = c("focal.base", "measure")) %>% 
                spread(key = measure, value = value) %>% 
                filter(base.call == focal.base) %>% 
                group_by(focal.base) %>% 
                summarize(avg.percsignal = mean(perc),
                    avg.areasignal = mean(area))
            
            # getting the critical value for each base
            crit.vals <- c(
                a = qZAGA(p = p.val.cutoff, mu = null.m.params$a$mu,
                    sigma = null.m.params$a$sigma,
                    nu = null.m.params$a$nu,
                    lower.tail = FALSE),
                c = qZAGA(p = p.val.cutoff, mu = null.m.params$c$mu,
                    sigma = null.m.params$c$sigma,
                    nu = null.m.params$c$nu,
                    lower.tail = FALSE),
                g = qZAGA(p = p.val.cutoff, mu = null.m.params$g$mu,
                    sigma = null.m.params$g$sigma,
                    nu = null.m.params$g$nu,
                    lower.tail = FALSE),
                t = qZAGA(p = p.val.cutoff, mu = null.m.params$t$mu,
                    sigma = null.m.params$t$sigma,
                    nu = null.m.params$t$nu,
                    lower.tail = FALSE)
            )
            # getting fillibens
            fil <- lapply(null.m.params, FUN = function(x){x$fillibens})
            filvec <- c(a = fil$a, c = fil$c, g = fil$g, t = fil$t)
            
            # getting mu -- a measure of dispersion
            mul <- lapply(null.m.params, FUN = function(x){x$mu})
            mulvec <- c(a = mul$a, c = mul$c, g = mul$g, t = mul$t)
            
            return(data.frame(avg.base, crit.perc.area = crit.vals, mu = mulvec, fillibens = filvec))
        })
        base.infoReactiveScramble[[i]] <- reactive({
            # this reactive is collecting information on each base:
            # - average percent signal
            # critical value for each base for the cutoff of significance
            
            sangs.filt <- sangs.filtReactive()[[input$nfiles+1]]
            null.m.params <- nullparams.ReactiveScramble()[[i]]
            p.val.cutoff <- p.val.Reactive()
            
            # finding the average percent signal for each base
            avg.base <- sangs.filt %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc) %>%
                separate(col = focal.base, into = c("focal.base", "measure")) %>% 
                spread(key = measure, value = value) %>% 
                filter(base.call == focal.base) %>% 
                group_by(focal.base) %>% 
                summarize(avg.percsignal = mean(perc),
                    avg.areasignal = mean(area))
            
            # getting the critical value for each base
            crit.vals <- c(
                a = qZAGA(p = p.val.cutoff, mu = null.m.params$a$mu,
                    sigma = null.m.params$a$sigma,
                    nu = null.m.params$a$nu,
                    lower.tail = FALSE),
                c = qZAGA(p = p.val.cutoff, mu = null.m.params$c$mu,
                    sigma = null.m.params$c$sigma,
                    nu = null.m.params$c$nu,
                    lower.tail = FALSE),
                g = qZAGA(p = p.val.cutoff, mu = null.m.params$g$mu,
                    sigma = null.m.params$g$sigma,
                    nu = null.m.params$g$nu,
                    lower.tail = FALSE),
                t = qZAGA(p = p.val.cutoff, mu = null.m.params$t$mu,
                    sigma = null.m.params$t$sigma,
                    nu = null.m.params$t$nu,
                    lower.tail = FALSE)
            )
            # getting fillibens
            fil <- lapply(null.m.params, FUN = function(x){x$fillibens})
            filvec <- c(a = fil$a, c = fil$c, g = fil$g, t = fil$t)
            
            # getting mu -- a measure of dispersion
            mul <- lapply(null.m.params, FUN = function(x){x$mu})
            mulvec <- c(a = mul$a, c = mul$c, g = mul$g, t = mul$t)
            
            return(data.frame(avg.base, crit.perc.area = crit.vals, mu = mulvec, fillibens = filvec))
        })
        
        output[[paste0('baseinfo.table',i)]] <- renderTable({
            base.info <- base.infoReactive[[i]]()
            
            temp <- base.info
            names(temp) <- c("Base", "Average percent signal",
                "Average peak area",  "Critical percent value",
                "model mu",  "Fillibens correlation")
            row.names(temp) <- NULL
            
            return(temp)
        })
        output[[paste0('baseinfo.tableScramble',i)]] <- renderTable({
            base.info <- base.infoReactiveScramble[[i]]()
            
            temp <- base.info
            names(temp) <- c("Base", "Average percent signal",
                "Average peak area",  "Critical percent value",
                "model mu",  "Fillibens correlation")
            row.names(temp) <- NULL
            
            return(temp)
        })
        
        ############################################################################
        # MAKE PREDICTED EDITING PLOTS
        ############################################################################
        
        output[[paste0('chromatogram_two',i)]] <- renderPlot({
            guide.coord <- guide.coordReactive()[[i]]
            input.basecalls <- input.basecallsReactive()[[i]]
            chromatogram(obj = input.basecalls,
                showcalls = "none",
                showhets = FALSE,
                trim5 = (guide.coord$start-1), trim3 = length(input.basecalls@primarySeq) - guide.coord$end,
                width = (guide.coord$end - guide.coord$start + 1))
        })
        output[[paste0('chromatogram_twoScramble',i)]] <- renderPlot({
            guide.coord <- guide.coordReactiveScramble()[[i]]
            input.basecallsScramble <- input.basecallsReactive()[[input$nfiles+1]]
            chromatogram(obj = input.basecallsScramble,
                showcalls = "none",
                showhets = FALSE,
                trim5 = (guide.coord$start-1), trim3 = length(input.basecallsScramble@primarySeq) - guide.coord$end,
                width = (guide.coord$end - guide.coord$start + 1))
        })
        
        output[[paste0('editing.table.plot',i)]] <- renderPlot({
            editing.df <- editing.Reactive()[[i]]
            sangs.filt <- sangs.filtReactive()[[i]]
            null.m.params <- nullparams.Reactive()[[i]]
            p.val.cutoff <- p.val.Reactive()
            
            #### Repeat code for getting avg.base from base.infoReactive
            # finding the average percent signal for each base
            avg.base <- sangs.filt %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc) %>%
                separate(col = focal.base, into = c("focal.base", "measure")) %>% 
                spread(key = measure, value = value) %>% 
                filter(base.call == focal.base) %>% 
                group_by(focal.base) %>% 
                summarize(avg.percsignal = mean(perc),
                    avg.areasignal = mean(area))
            
            # finding the model mu
            mul <- lapply(null.m.params, FUN = function(x){x$mu})
            mulvec <- c(a = mul$a, c = mul$c, g = mul$g, t = mul$t)
            ####
            
            ### Reshape data
            
            edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
                separate(col = focal.base, into = c("focal.base", "measure"))
            
            edit.spread <- edit.long %>% 
                spread(key = measure, value = value) 
            
            color.cutoff = min(avg.base$avg.percsignal - mulvec)
            edit.color <- edit.spread %>% 
                mutate(adj.perc = {ifelse(perc >= color.cutoff,
                    100,
                    perc)
                } %>% as.numeric) %>%
                filter(pval < p.val.cutoff)
            
            
            
            #### make editing_table
            if(any(edit.color$adj.perc != 100)){
                edit.spread %>%
                    ggplot(aes(x = as.factor(index), y = focal.base)) + 
                    geom_tile(data = edit.color, aes(fill = adj.perc)) + 
                    geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +   
                    guides(fill = FALSE) + 
                    scale_fill_continuous(low = "#f7a8a8", high = "#9acdee") + 
                    scale_x_discrete(position = "top", labels = editing.df$guide.seq) + 
                    labs(x = NULL, y = NULL) + 
                    theme(axis.ticks = element_blank(),
                        axis.text=element_text(size=16),
                        plot.title = element_text(hjust = 0, size = 16),
                        plot.margin=unit(c(0,0,0,2), "cm"), #c(top, bottom, left, right)
                        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                        plot.background = element_rect(fill = "transparent",colour = NA)
                    ) +
                    coord_fixed(1)} 
            else
            {edit.spread %>%
                    ggplot(aes(x = as.factor(index), y = focal.base)) + 
                    geom_tile(data = edit.color, fill = "#9acdee") + 
                    geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +   
                    guides(fill = FALSE) + 
                    scale_x_discrete(position = "top", labels = editing.df$guide.seq) + 
                    labs(x = NULL, y = NULL) + 
                    theme(axis.ticks = element_blank(),
                        axis.text=element_text(size=16),
                        plot.title = element_text(hjust = 0, size = 16),
                        plot.margin=unit(c(0,0,0,2), "cm"), #c(top, bottom, left, right)
                        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                        plot.background = element_rect(fill = "transparent",colour = NA)
                    ) +
                    coord_fixed(1)}
        })
        output[[paste0('editing.table.plotScramble',i)]] <- renderPlot({
            editing.df <- editing.ReactiveScramble()[[i]]
            sangs.filt <- sangs.filtReactive()[[input$nfiles+1]]
            null.m.params <- nullparams.ReactiveScramble()[[i]]
            p.val.cutoff <- p.val.Reactive()
            
            #### Repeat code for getting avg.base from base.infoReactive
            # finding the average percent signal for each base
            avg.base <- sangs.filt %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc) %>%
                separate(col = focal.base, into = c("focal.base", "measure")) %>% 
                spread(key = measure, value = value) %>% 
                filter(base.call == focal.base) %>% 
                group_by(focal.base) %>% 
                summarize(avg.percsignal = mean(perc),
                    avg.areasignal = mean(area))
            
            # finding the model mu
            mul <- lapply(null.m.params, FUN = function(x){x$mu})
            mulvec <- c(a = mul$a, c = mul$c, g = mul$g, t = mul$t)
            ####
            
            ### Reshape data
            
            edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
                separate(col = focal.base, into = c("focal.base", "measure"))
            
            edit.spread <- edit.long %>% 
                spread(key = measure, value = value) 
            
            color.cutoff = min(avg.base$avg.percsignal - mulvec)
            edit.color <- edit.spread %>% 
                mutate(adj.perc = {ifelse(perc >= color.cutoff,
                    100,
                    perc)
                } %>% as.numeric) %>%
                filter(pval < p.val.cutoff)
            
            
            
            #### make editing_table
            if(any(edit.color$adj.perc != 100)){
                edit.spread %>%
                    ggplot(aes(x = as.factor(index), y = focal.base)) + 
                    geom_tile(data = edit.color, aes(fill = adj.perc)) + 
                    geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +   
                    guides(fill = FALSE) + 
                    scale_fill_continuous(low = "#f7a8a8", high = "#9acdee") + 
                    scale_x_discrete(position = "top", labels = editing.df$guide.seq) + 
                    labs(x = NULL, y = NULL) + 
                    theme(axis.ticks = element_blank(),
                        axis.text=element_text(size=16),
                        plot.title = element_text(hjust = 0, size = 16),
                        plot.margin=unit(c(0,0,0,2), "cm"), #c(top, bottom, left, right)
                        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                        plot.background = element_rect(fill = "transparent",colour = NA)
                    ) +
                    coord_fixed(1)} 
            else
            {edit.spread %>%
                    ggplot(aes(x = as.factor(index), y = focal.base)) + 
                    geom_tile(data = edit.color, fill = "#9acdee") + 
                    geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +   
                    guides(fill = FALSE) + 
                    scale_x_discrete(position = "top", labels = editing.df$guide.seq) + 
                    labs(x = NULL, y = NULL) + 
                    theme(axis.ticks = element_blank(),
                        axis.text=element_text(size=16),
                        plot.title = element_text(hjust = 0, size = 16),
                        plot.margin=unit(c(0,0,0,2), "cm"), #c(top, bottom, left, right)
                        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                        plot.background = element_rect(fill = "transparent",colour = NA)
                    ) +
                    coord_fixed(1)}
        })
        
        ############################################################################
        # DOWNLOAD FULL PLOT
        ############################################################################
        
        editPlot[[i]] <- reactive({
            editing.df <- editing.Reactive()[[i]]
            sangs.filt <- sangs.filtReactive()[[i]]
            null.m.params <- nullparams.Reactive()[[i]]
            p.val.cutoff <- p.val.Reactive()
            
            #### Repeat code for getting avg.base from base.infoReactive
            # finding the average percent signal for each base
            avg.base <- sangs.filt %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc) %>%
                separate(col = focal.base, into = c("focal.base", "measure")) %>% 
                spread(key = measure, value = value) %>% 
                filter(base.call == focal.base) %>% 
                group_by(focal.base) %>% 
                summarize(avg.percsignal = mean(perc),
                    avg.areasignal = mean(area))
            
            # finding the model mu
            mul <- lapply(null.m.params, FUN = function(x){x$mu})
            mulvec <- c(a = mul$a, c = mul$c, g = mul$g, t = mul$t)
            ####
            
            ### Reshape data
            
            edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
                separate(col = focal.base, into = c("focal.base", "measure"))
            
            edit.spread <- edit.long %>% 
                spread(key = measure, value = value) 
            
            color.cutoff = min(avg.base$avg.percsignal - mulvec)
            edit.color <- edit.spread %>% 
                mutate(adj.perc = {ifelse(perc >= color.cutoff,
                    100,
                    perc)
                } %>% as.numeric) %>%
                filter(pval < p.val.cutoff)
            
            
            
            #### make editing_table
            if(any(edit.color$adj.perc != 100)){
                edit.spread %>%
                    ggplot(aes(x = as.factor(index), y = focal.base)) + 
                    geom_tile(data = edit.color, aes(fill = adj.perc)) + 
                    geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +   
                    guides(fill = FALSE) + 
                    scale_fill_continuous(low = "#f7a8a8", high = "#9acdee") + 
                    scale_x_discrete(position = "top", labels = editing.df$guide.seq) + 
                    labs(x = NULL, y = NULL) + 
                    theme(axis.ticks = element_blank(),
                        axis.text=element_text(size=16),
                        plot.title = element_text(hjust = 0, size = 16),
                        plot.margin=unit(c(0,0,0,2), "cm"), #c(top, bottom, left, right)
                        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                        plot.background = element_rect(fill = "transparent",colour = NA)
                    ) +
                    coord_fixed(1)} 
            else
            {edit.spread %>%
                    ggplot(aes(x = as.factor(index), y = focal.base)) + 
                    geom_tile(data = edit.color, fill = "#9acdee") + 
                    geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +   
                    guides(fill = FALSE) + 
                    scale_x_discrete(position = "top", labels = editing.df$guide.seq) + 
                    labs(x = NULL, y = NULL) + 
                    theme(axis.ticks = element_blank(),
                        axis.text=element_text(size=16),
                        plot.title = element_text(hjust = 0, size = 16),
                        plot.margin=unit(c(0,0,0,2), "cm"), #c(top, bottom, left, right)
                        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                        plot.background = element_rect(fill = "transparent",colour = NA)
                    ) +
                    coord_fixed(1)}
        })
        editPlotScramble[[i]] <- reactive({
            editing.df <- editing.ReactiveScramble()[[i]]
            sangs.filt <- sangs.filtReactive()[[input$nfiles+1]]
            null.m.params <- nullparams.ReactiveScramble()[[i]]
            p.val.cutoff <- p.val.Reactive()
            
            #### Repeat code for getting avg.base from base.infoReactive
            # finding the average percent signal for each base
            avg.base <- sangs.filt %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc) %>%
                separate(col = focal.base, into = c("focal.base", "measure")) %>% 
                spread(key = measure, value = value) %>% 
                filter(base.call == focal.base) %>% 
                group_by(focal.base) %>% 
                summarize(avg.percsignal = mean(perc),
                    avg.areasignal = mean(area))
            
            # finding the model mu
            mul <- lapply(null.m.params, FUN = function(x){x$mu})
            mulvec <- c(a = mul$a, c = mul$c, g = mul$g, t = mul$t)
            ####
            
            ### Reshape data
            
            edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
                separate(col = focal.base, into = c("focal.base", "measure"))
            
            edit.spread <- edit.long %>% 
                spread(key = measure, value = value) 
            
            color.cutoff = min(avg.base$avg.percsignal - mulvec)
            edit.color <- edit.spread %>% 
                mutate(adj.perc = {ifelse(perc >= color.cutoff,
                    100,
                    perc)
                } %>% as.numeric) %>%
                filter(pval < p.val.cutoff)
            
            
            
            #### make editing_table
            if(any(edit.color$adj.perc != 100)){
                edit.spread %>%
                    ggplot(aes(x = as.factor(index), y = focal.base)) + 
                    geom_tile(data = edit.color, aes(fill = adj.perc)) + 
                    geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +   
                    guides(fill = FALSE) + 
                    scale_fill_continuous(low = "#f7a8a8", high = "#9acdee") + 
                    scale_x_discrete(position = "top", labels = editing.df$guide.seq) + 
                    labs(x = NULL, y = NULL) + 
                    theme(axis.ticks = element_blank(),
                        axis.text=element_text(size=16),
                        plot.title = element_text(hjust = 0, size = 16),
                        plot.margin=unit(c(0,0,0,2), "cm"), #c(top, bottom, left, right)
                        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                        plot.background = element_rect(fill = "transparent",colour = NA)
                    ) +
                    coord_fixed(1)} 
            else
            {edit.spread %>%
                    ggplot(aes(x = as.factor(index), y = focal.base)) + 
                    geom_tile(data = edit.color, fill = "#9acdee") + 
                    geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +   
                    guides(fill = FALSE) + 
                    scale_x_discrete(position = "top", labels = editing.df$guide.seq) + 
                    labs(x = NULL, y = NULL) + 
                    theme(axis.ticks = element_blank(),
                        axis.text=element_text(size=16),
                        plot.title = element_text(hjust = 0, size = 16),
                        plot.margin=unit(c(0,0,0,2), "cm"), #c(top, bottom, left, right)
                        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                        plot.background = element_rect(fill = "transparent",colour = NA)
                    ) +
                    coord_fixed(1)}
        })
        
        output[[paste0('downloadMulti',i)]] <-  downloadHandler(
            # download data
            filename = function() {
                n <- filename_out()[[i]]
                return(paste0(n, '.pdf'))},
            content = function(file) {
                guide.coord <- guide.coordReactive()[[i]]
                input.basecalls <- input.basecallsReactive()[[i]]
                temp_chrom <- file.path(tempdir(), paste0("chrom",i,".pdf"))
                chromatogram(obj = input.basecalls,
                    showcalls = "none",
                    showhets = FALSE,
                    trim5 = (guide.coord$start-1), trim3 = length(input.basecalls@primarySeq) - guide.coord$end,
                    width = (guide.coord$end - guide.coord$start + 1), filename=temp_chrom)
                
                guide.coordScramble <- guide.coordReactiveScramble()[[i]]
                input.basecallsScramble <- input.basecallsReactive()[[input$nfiles+1]]
                temp_chromScramble  <- file.path(tempdir(), paste0("chromScramble",i,".pdf"))
                chromatogram(obj = input.basecallsScramble,
                    showcalls = "none",
                    showhets = FALSE,
                    trim5 = (guide.coordScramble$start-1), trim3 = length(input.basecallsScramble@primarySeq) - guide.coordScramble$end,
                    width = (guide.coordScramble$end - guide.coordScramble$start + 1), filename=temp_chromScramble)
                
                p1 <- ggdraw() + draw_image(magick::image_read_pdf(temp_chrom, density = 300),hjust =0.01)
                p2 <- editPlot[[i]]()
                p3 <- ggdraw() + draw_image(magick::image_read_pdf(temp_chromScramble, density = 300),hjust=0.01)
                p4 <- editPlotScramble[[i]]()
                p5 <- plot_grid(p1,p2,p3,p4,ncol=1,axis="l", align="v", scale=c(1.15,1,1.15,1), labels=c(paste0(
                    "Guide: ", input[[paste0('guide', i)]]), "E", "Scramble:", "F"), label_x=c(-0.2,-1,0,-1))
                pdf(file)
                print(p5)
                dev.off()
            }
        )
        
        ############################################################################
        # MAKE QUAD PLOT
        ############################################################################
        output[[paste0('editing.quad.plot',i)]] <- renderPlot({
            editing.df <- editing.Reactive()[[i]]
            null.m.params <- nullparams.Reactive()[[i]]
            p.val.cutoff <- p.val.Reactive()
            
            edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
                separate(col = focal.base, into = c("focal.base", "measure"))
            
            p.a <- makeEditingBarPlot(edit.long = edit.long, null.m.params = null.m.params$a,
                base = "A", pval = p.val.cutoff, editing.df)
            p.c <- makeEditingBarPlot(edit.long, null.m.params$c,
                base = "C", pval = p.val.cutoff, editing.df)
            p.g <- makeEditingBarPlot(edit.long, null.m.params$g,
                base = "G", pval = p.val.cutoff, editing.df)
            p.t <- makeEditingBarPlot(edit.long, null.m.params$t,
                base = "T", pval = p.val.cutoff, editing.df)
            
            grid.arrange(p.a, p.c, p.g, p.t)
        })
        output[[paste0('editing.quad.plotScramble',i)]] <- renderPlot({
            editing.df <- editing.ReactiveScramble()[[i]]
            null.m.params <- nullparams.ReactiveScramble()[[i]]
            p.val.cutoff <- p.val.Reactive()
            
            edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
                separate(col = focal.base, into = c("focal.base", "measure"))
            
            p.a <- makeEditingBarPlot(edit.long = edit.long, null.m.params = null.m.params$a,
                base = "A", pval = p.val.cutoff, editing.df)
            p.c <- makeEditingBarPlot(edit.long, null.m.params$c,
                base = "C", pval = p.val.cutoff, editing.df)
            p.g <- makeEditingBarPlot(edit.long, null.m.params$g,
                base = "G", pval = p.val.cutoff, editing.df)
            p.t <- makeEditingBarPlot(edit.long, null.m.params$t,
                base = "T", pval = p.val.cutoff, editing.df)
            
            grid.arrange(p.a, p.c, p.g, p.t)
        })
        
        ############################################################################
        # MAKE TABLE OF EDITING RESULTS (REPORT)
        ############################################################################
        base.editingtabledata[[i]] <- reactive({ # maybe need to use renderMarkdown?
            editing.df <- editing.Reactive()[[i]]
            p.val.cutoff <- p.val.Reactive()
            
            edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
                separate(col = focal.base, into = c("focal.base", "measure"))
            
            edit.spread <- edit.long %>% spread(key = measure, value = value)
            
            
            editingtable <- edit.spread %>% arrange(guide.position) %>% select(index,guide.position, guide.seq, base.call, focal.base, perc, pval) %>% 
                mutate(index = as.character(index), 
                    guide.position = as.character(guide.position),
                    perc = format(perc, digits = 2),
                    signficant = ifelse(pval < p.val.cutoff, yes = "*", no = ""))
            
            names(editingtable) <-  c("Sanger_position", "Guide_position","Guide_sequence", "Sanger_base_call",
                "Focal_base", "Focal_base_peak_area", "p_value", "Significance")
            
            return(editingtable)
        })
        base.editingtabledataScramble[[i]] <- reactive({ # maybe need to use renderMarkdown?
            editing.df <- editing.ReactiveScramble()[[i]]
            p.val.cutoff <- p.val.Reactive()
            
            edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
                separate(col = focal.base, into = c("focal.base", "measure"))
            
            edit.spread <- edit.long %>% spread(key = measure, value = value)
            
            
            editingtable <- edit.spread %>% arrange(guide.position) %>% select(index,guide.position, guide.seq, base.call, focal.base, perc, pval) %>% 
                mutate(index = as.character(index), 
                    guide.position = as.character(guide.position),
                    perc = format(perc, digits = 2),
                    signficant = ifelse(pval < p.val.cutoff, yes = "*", no = ""))
            
            names(editingtable) <-  c("Sanger_position", "Guide_position","Guide_sequence", "Sanger_base_call",
                "Focal_base", "Focal_base_peak_area", "p_value", "Significance")
            
            return(editingtable)
        })
        
        output[[paste0('editingtable',i)]] <- renderTable(digits = -3, expr = base.editingtabledata[[i]]())
        output[[paste0('editingtableScramble',i)]] <- renderTable(digits = -3, expr = base.editingtabledataScramble[[i]]())
        
        output[[paste0('modelfits',i)]] <- renderTable({
            null.m.params <- nullparams.Reactive()[[i]]
            
            temp <- data.frame(`Base` = c("A", "C", "G", "T"),
                `Fillibens Correlation` = c(null.m.params$a$fillibens, 
                    null.m.params$g$fillibens, 
                    null.m.params$c$fillibens, 
                    null.m.params$t$fillibens))
            
            names(temp) <- c("Base", "Fillibens Correlation")
            return(temp)
        })
        output[[paste0('modelfitsScramble',i)]] <- renderTable({
            null.m.params <- nullparams.ReactiveScramble()[[i]]
            
            temp <- data.frame(`Base` = c("A", "C", "G", "T"),
                `Fillibens Correlation` = c(null.m.params$a$fillibens, 
                    null.m.params$g$fillibens, 
                    null.m.params$c$fillibens, 
                    null.m.params$t$fillibens))
            
            names(temp) <- c("Base", "Fillibens Correlation")
            return(temp)
        })
        
        ############################################################################
        # REPORT DOWNLOAD
        ############################################################################
        
        base.preocessdData[[i]] <- reactive({
            data <- base.editingtabledata[[i]]()
            dataScramble <- base.editingtabledataScramble[[i]]()
            rev.guide <- eval(parse(text=paste0('input$guide.rev', i)))
            p.val.cutoff <- p.val.Reactive()
            if (input$editorType=='ABE') {
                edited_base <- 'A'
                edited_base_rev <- 'T'
                focal_base <- 'G'
                focal_base_rev <- 'C'
            } else if  (input$editorType=='CBE') {
                edited_base <- 'C'
                edited_base_rev <- 'G'
                focal_base <- 'T'
                focal_base_rev <- 'A'
            }
            if(!rev.guide) {
                filtered_data <- data %>% filter(Guide_sequence == edited_base) %>% filter(Focal_base == focal_base) %>% select(
                    c("Guide_position", "Focal_base_peak_area"))
                filtered_dataScramble <- dataScramble %>% filter(Guide_sequence == edited_base) %>% filter(
                    Focal_base == focal_base) %>% select(c("Guide_position", "Focal_base_peak_area"))
            } else {
                filtered_data <- data %>% filter(Guide_sequence == edited_base_rev) %>% filter(Focal_base == focal_base_rev) %>% select(
                    c("Guide_position", "Focal_base_peak_area"))
                filtered_dataScramble <- dataScramble %>% filter(Guide_sequence == edited_base_rev) %>% filter(
                    Focal_base == focal_base_rev) %>% select(c("Guide_position", "Focal_base_peak_area"))
            }
            colnames(filtered_dataScramble) <- c("Guide_position", "Focal_base_peak_area_scramble")
            filtered_data <- filtered_data %>% full_join(filtered_dataScramble)
            filtered_data$Difference <- as.double(filtered_data[["Focal_base_peak_area"]]) - as.double(
                filtered_data[["Focal_base_peak_area_scramble"]])
            filtered_data$Difference <- sapply(filtered_data$Difference, function(x) max(0, x))
            
            if((!rev.guide & input$orientation==3) | (rev.guide & input$orientation==5)) {
                filtered_data$A <- nchar(eval(parse(text=paste0('input$guide', i))))-as.integer(filtered_data$Guide_position)+1
            } else {
                filtered_data$A <- as.integer(filtered_data$Guide_position)
            }
            
            filtered_data$A <- sapply(filtered_data$A, function(x) paste0(edited_base,x))
            filtered_data <- filtered_data %>% select(c("A", "Focal_base_peak_area", "Focal_base_peak_area_scramble", "Difference"))
            colnames(filtered_data) <- c(paste0(edited_base, "# (from ", input$orientation, "')"), "Guide", "Scramble", "Difference")
            if(((!rev.guide & input$orientation==3) | (rev.guide & input$orientation==5)) & dim(filtered_data)[1]>1) {
                rev_data_frame <- apply(filtered_data, 2, rev)
                return(tibble(as.data.frame(rev_data_frame)))
            } else {
                return(filtered_data)
            }
        })
        env[['base.preocessdDataEnv']][[i]] <- reactive({
            data <- base.editingtabledata[[i]]()
            dataScramble <- base.editingtabledataScramble[[i]]()
            rev.guide <- eval(parse(text=paste0('input$guide.rev', i)))
            p.val.cutoff <- p.val.Reactive()
            if (input$editorType=='ABE') {
                edited_base <- 'A'
                edited_base_rev <- 'T'
                focal_base <- 'G'
                focal_base_rev <- 'C'
            } else if  (input$editorType=='CBE') {
                edited_base <- 'C'
                edited_base_rev <- 'G'
                focal_base <- 'T'
                focal_base_rev <- 'A'
            }
            if(!rev.guide) {
                filtered_data <- data %>% filter(Guide_sequence == edited_base) %>% filter(Focal_base == focal_base) %>% select(
                    c("Guide_position", "Focal_base_peak_area"))
                filtered_dataScramble <- dataScramble %>% filter(Guide_sequence == edited_base) %>% filter(
                    Focal_base == focal_base) %>% select(c("Guide_position", "Focal_base_peak_area"))
            } else {
                filtered_data <- data %>% filter(Guide_sequence == edited_base_rev) %>% filter(Focal_base == focal_base_rev) %>% select(
                    c("Guide_position", "Focal_base_peak_area"))
                filtered_dataScramble <- dataScramble %>% filter(Guide_sequence == edited_base_rev) %>% filter(
                    Focal_base == focal_base_rev) %>% select(c("Guide_position", "Focal_base_peak_area"))
            }
            colnames(filtered_dataScramble) <- c("Guide_position", "Focal_base_peak_area_scramble")
            filtered_data <- filtered_data %>% full_join(filtered_dataScramble)
            filtered_data$Difference <- as.double(filtered_data[["Focal_base_peak_area"]]) - as.double(
                filtered_data[["Focal_base_peak_area_scramble"]])
            filtered_data$Difference <- sapply(filtered_data$Difference, function(x) max(0, x))
            
            if((!rev.guide & input$orientation==3) | (rev.guide & input$orientation==5)) {
                filtered_data$A <- nchar(eval(parse(text=paste0('input$guide', i))))-as.integer(filtered_data$Guide_position)+1
            } else {
                filtered_data$A <- as.integer(filtered_data$Guide_position)
            }
            
            filtered_data$A <- sapply(filtered_data$A, function(x) paste0(edited_base,x))
            filtered_data <- filtered_data %>% select(c("A", "Focal_base_peak_area", "Focal_base_peak_area_scramble", "Difference"))
            colnames(filtered_data) <- c(paste0(edited_base, "# (from ", input$orientation, "')"), "Guide", "Scramble", "Difference")
            if(((!rev.guide & input$orientation==3) | (rev.guide & input$orientation==5)) & dim(filtered_data)[1]>1) {
                rev_data_frame <- apply(filtered_data, 2, rev)
                return(tibble(as.data.frame(rev_data_frame)))
            } else {
                return(filtered_data)
            }
        })
        
        output[[paste0('preocessdData.table',i)]] <- renderTable({
            temp <- base.preocessdData[[i]]()
            row.names(temp) <- NULL
            return(temp)
        })
        
        output[[paste0('downloadData',i)]] <-  downloadHandler(
            # download data
            filename = function() {
                n <- filename_out()[[i]]
                return(paste0(n, '.csv'))},
            content = function(file) {
                # save data
                write.table(base.preocessdData[[i]](), file = file, quote = FALSE, sep='\t', row.names = FALSE)
            }
        )
        
        output[[paste0('samplename1',i)]] <- reactive(paste0('Sample: ', filename_out()[[i]]))
        output[[paste0('samplename2',i)]] <- reactive(paste0('Sample: ', filename_out()[[i]]))
        output[[paste0('samplename3',i)]] <- reactive(paste0('Sample: ', filename_out()[[i]]))
        
    })
    
    # make excle file with all data
    output$downloadAllData <- downloadHandler(
        # download data
        filename = 'bacth_data.xlsx',
        content = function(file) {
            # save data
            wb <- createWorkbook()
            sheetnames <- filename_out()
            data <- lapply(1:input$nfiles, function (i) base.preocessdDataEnv[[i]]())
            sheets <- lapply(sheetnames, createSheet, wb = wb)
            void <- Map(addDataFrame, data, sheets)
            saveWorkbook(wb, file = file)
        }
    )
    })
    
    ############################################################################
    # POPULATE OUTPUT TABS
    ############################################################################
    
    output$outputQC <- renderUI({
        req(input$nfiles)
        output_list <- list(
            h2("Data QC"),
            p('Data QC for each sample and control pair.')
        )
        output_list <- lapply(1:input$nfiles, function(i) {
            output_list_i <- list(
                h3(textOutput(paste0("samplename1", i))),
                h4("Total peak area before filtering"),
                plotOutput(outputId = paste0("prefilter.totalarea",i)),
                h4("Scramble:"),
                plotOutput(outputId = paste0("prefilter.totalareaScramble",i)),
                p("This plot shows the total peak area at each positon before QC filtering.
                The peak are for each base (ACGT) was summed to produce the total
                peak area. The shaded region denotes the guide RNA"),
                
                h4("Data QA: Signal and noise plot"),
                plotlyOutput(outputId =  paste0("postfilter.signal.noise",i)),
                h4("Scramble:"),
                plotlyOutput(outputId = paste0("postfilter.signal.noiseScramble",i)),
                textOutput(outputId =  paste0("trimmedrange", i)),
                br(),
                p(strong("Signal and noise plot:"), " This plot shows the peak area of the filtered dataset. This can tell
                you if there are particular regions in your sequencing that have high
                noise, and that it might affect the sensitivity of the predicted
                editing that we can make. The shaded region denotes the guide RNA region.
                This plot shows the data postfiltering that will be used to generate the null
                distribtution for our prediction method. The first 20 bases are removed as they
                are generally poor quality. Then we set a filtering cutoff based on the mean of
                the toal peak area, we remove positions where the total peak area is less than one
                tenth of the mean. For example, if the mean total peak area is 1000, we remove positions
                that have less than 100 total peak area from our dataset. This mainly serves to trim the
                low quality bases at the tail end of the sequencing run."),
                
                h4("Data QA: Percent noise peak area"),
                plotOutput(outputId =  paste0("postfilter.noise.perc",i)),
                h4("Scramble:"),
                plotOutput(outputId =  paste0("postfilter.noise.percScramble",i)),
                p("This plot shows the amount of noise at a base position as the percent of the total peak area.
                The shaded region is where the gRNA matches, which might have a peak due to base editing."),
                br()
            )
            output_list <- c(output_list, output_list_i)
        })
        do.call(tagList, output_list)
    })
    
    output$outputEditing <- renderUI({
        req(input$nfiles)
        output_list <- list(
            h2("Predicted Editing"),
            p('Predicted editing for each sample and control pair.')
        )
        output_list <- lapply(1:input$nfiles, function(i) {
            output_list_i <- list(
                h3(textOutput(paste0("samplename2", i))),
                h4("Predicted Editing"),
                p("Guide:"),
                textOutput(outputId = paste0("guideSeq",i)),
                fluidRow(
                    verticalLayout(
                        plotOutput(outputId =  paste0("chromatogram_two",i)),
                        plotOutput(outputId =  paste0("editing.table.plot",i), width = "93.5%"))),
                h4("Scramble"),
                fluidRow(
                    verticalLayout(
                        plotOutput(outputId = paste0("chromatogram_twoScramble",i)),
                        plotOutput(outputId = paste0("editing.table.plotScramble",i), width = "93.5%"))),
                downloadButton( paste0('downloadMulti',i)),
                p("The top plot is the chromatogram of the protospacer. Highlighted peaks indicate double peaks were detected, however a peak may still be significant even if not highlighted. \nThe bottom plot shows the percent area of the signal for each base (ACGT) at each position along
            the guide. Bases that are significantly different from the noise are colored in, color coded
            relative to their percent area. Most positions of the guide only have one base colored in, as there
            this is due to only one peak being present at that position. In the case of editing, there are more
            than one base colored in at a position."),
                h4("Base info"),
                p("You can use the information in the following table to better understand the values given for the percent
            area in the 'Table plot'. The average percent signal shows you the average percent area (area of signal / all
            other base areas), and (100 - average percent signal) shows you how much noise is present during the reading
            of that base. The amount of noise present in each base is also reflected in the 'Critical percent value',
            which is the critical value for the percent area for that base -- any percent area value above that would
            be called significantly different from the noise, and therefore editing has occurred. This critical value is
            calculated based on the P-value cutoff specified (default is 0.01). The model mu is the mu term for the zero
            adjusted gamma distribution to model the noise -- a higher value of mu means that the predicted editing value
            is less accurate. Filliben's correlation is how well the noise is modelled by a zero adjusted gamma
            distribution. Lower values of Filliben's correlation means less confidence in EditR for predicting editing,
            your values should be above 0.90, if not, this means that your sequencing file likely has problems
            with data quality."),
                tableOutput(outputId =  paste0("baseinfo.table",i)),
                h4("Scramble:"),
                tableOutput(outputId = paste0("baseinfo.tableScramble",i)),
                h4("Quad plot"),
                plotOutput(outputId =  paste0("editing.quad.plot",i)),
                h4("Scramble:"),
                plotOutput(outputId = paste0("editing.quad.plotScramble",i)),
                p("This plot shows the percent area of the background bases (the bases that
            were not the guide squence), the line denotes the crtical value cutoff based on the P-value
            cutoff specified (default is 0.010.")
            )
            output_list <- c(output_list, output_list_i)
        })
        do.call(tagList, output_list)
    })
    
    output$downloadData <- renderUI({
        req(input$nfiles)
        output_list <- list(
            h2("Download Data"),
            p('Download data for each sample and control pair.')
        )
        output_list <- lapply(1:input$nfiles, function(i) {
            output_list_i <- list(
                h3(textOutput(paste0("samplename3", i))),
                h4("Base editing data"),
                tableOutput(outputId = paste0("preocessdData.table",i)),
                p("Download base editing data"),
                downloadButton(paste0('downloadData',i))
            )
            output_list <- c(output_list, output_list_i)
        })
        do.call(tagList, output_list)
    })
    

    
    })
