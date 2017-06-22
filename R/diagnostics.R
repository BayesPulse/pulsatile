# 
# 
# mcmc_trace <- function() {}
# mcmc_posteriors <- function() {}
# mcmc_locations <- function() {}
# 
# 
# 
#   temp <- 
#     common %>%
#       gather(key = key, value = value, num.pulses:sd.widths) %>%
#       filter(dataset %in% dataset.nums) %>%
#       group_by(dataset) %>%
#       do( 
#         # Trace plots
#         trace.figs = 
#         {
#           trace.fig <- 
#             ggplot(., aes(x = iteration, y = value)) +
#               geom_path(size = 0.10) +
#               facet_wrap( ~ key, ncol = 2, nrow = 4, scales = "free") +
#               ggtitle(paste("Dataset", ifelse(is.null(.$orig.dataset), 
#                                               unique(.$dataset),
#                                               unique(.$orig.dataset))))
#           #print(trace.fig)
#         },
#         # Posterior densities
#         post.figs = 
#         {
#           post.fig <- 
#             ggplot(., aes(x = value)) +
#               geom_histogram(aes(y = ..density..), size = 0.15) +
#               #geom_density(alpha=.2, fill="#FF6666") +
#               facet_wrap( ~ key, ncol = 2, nrow = 4, scales = "free") +
#               ggtitle(paste("Dataset", ifelse(is.null(.$orig.dataset), 
#                                               unique(.$dataset),
#                                               unique(.$orig.dataset))))
#           #suppressMessages(print(post.fig))
#         }
#       )
# 
#   location.figs.lst <- 
#     pulse %>%
#       filter(dataset %in% dataset.nums) %>%
#       group_by(dataset) %>%
#       do(
#         # Location histograms
#         location.figs = 
#         {
#           location.fig <- 
#             ggplot(., aes(x = location)) +
#               geom_histogram(binwidth = 5) +
#               theme(panel.grid.minor = element_line(colour="lightgrey", size=0.5)) + 
#               theme(panel.grid.major = element_line(colour="lightgrey", size=0.5)) + 
#               scale_x_continuous(breaks = seq(-50, max.time+50, 50),
#                                  minor_breaks = seq(-50, max.time+50, 10),
#                                  limits = c(-50, max.time+50)) 
#         }
#       )
# 
#   sim.figs.lst <-
#     sim %>%
#     filter(dataset %in% dataset.nums) %>%
#     group_by(dataset) %>%
#     do(
#        sim.figs = 
#        {
#          sim.fig <-
#            ggplot(., aes(x = time, y = concentration)) +
#              geom_path() +
#              geom_point() + 
#              theme(panel.grid.minor = element_line(colour="lightgrey", size=0.5)) + 
#              theme(panel.grid.major = element_line(colour="lightgrey", size=0.5)) + 
#              scale_x_continuous(breaks = seq(-50, max.time+50, 50),
#                                 minor_breaks = seq(-50, max.time+50, 10),
#                                 limits = c(-50, max.time+50)) 
#        }
#     )
