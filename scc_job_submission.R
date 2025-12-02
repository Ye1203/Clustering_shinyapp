submit_scc_job <- function(seuratObj,
                           save_path,
                           project_name,
                           start_resolution,
                           step_resolution,
                           end_resolution,
                           runtime,
                           cores,
                           email) {
  saveRDS(seuratObj, file = file.path(save_path, "input_seurat_file.rds"))
  r_script_path <- file.path(save_path, "scc_job_script.R")
  # Generate R script file  
  r_code <- sprintf(
    '
library(Seurat)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(cowplot)
source("renv/activate.R")

output_file_path <- "%s"
email_address <- "%s"
start_resolution <- "%s"
step_resolution <- "%s"
end_resolution <- "%s"
send_email <- function(to, subject, body) {
  f <- tempfile(fileext = ".txt")
  mail_content <- c(
    paste0("To: ", to),
    paste0("Subject: ", subject),
    "",
    body
  )
  writeLines(mail_content, f)
  cmd <- sprintf("/usr/sbin/sendmail -t < %%s", shQuote(f))  
  ret <- system(cmd)
  unlink(f)
  if (ret != 0) {
    warning("sendmail command failed with exit code ", ret)
  }
}
start_time <- Sys.time()

seurat_obj <- readRDS(file.path(output_file_path, "input_seurat_file.rds"))


  if (!is.na(email_address)) {
    body <- c(
      "Hi,",
      "",
      sprintf("The procedure has been completed successfully, the result is under folder %%s. Spending %%s %%s", output_file_path, round(as.numeric(elapsed_time),2), attr(elapsed_time, "units")),
      "",
      "Best,",
      "Bingtian"
    )
    send_email(to = email_address, subject = "Program Completed", body = body)
  }
  
  list(success = TRUE, elapsed = elapsed_time)
}, error = function(e) {
  err_msg <- conditionMessage(e)
  
  if (!is.na(email_address)) {
    body <- c(
      "Hi,",
      "",
      sprintf("There is an error (%%s). You can copy this email and send to btye@bu.edu.", err_msg),
      "",
      sprintf("output_file_path = %%s", output_file_path),
      "",
      "Best,",
      "Bingtian"
    )
    send_email(to = email_address, subject = "Program Completed", body = body)
  }
  line  <- sprintf("Spending %%s %%s", as.numeric(elapsed_time), attr(elapsed_time, "units"))
  writeLines(
  line,
  con = file.path(output_file_path, "time.txt")
)
  list(success = FALSE, error = err_msg)
})
',output_file_path,
  email_address,
  start_resolution,
  step_resolution,
  end_resolution
  )

write.Lines(r_code, con = r_script_path)

  # Generate qsub file
  qsub_file_path <- file.path(save_path, "launch_job.qsub")
  log_path <- file.path(save_path, "lamian_job.log")
  qsub_content <- paste0(
    "#!/bin/bash\n",
    "#$ -N clustering\n",
    "#$ -cwd\n",
    "#$ -j y\n",
    "#$ -o ", log_path, "\n",
    "#$ -pe omp ", cores, "\n",
    "#$ -l h_rt=", runtime, ":00:00\n",
    "#$ -V\n",
    "#$ -P ", project_name, "\n",
    "\n",
    "echo \"==========================================================\"\n",
    "Start_Time=$(date +\"%s\")\n",
    "echo \"Starting on       : $(date)\"\n",
    "echo \"Running on node   : $(hostname)\"\n",
    "echo \"Current job ID    : $JOB_ID\"\n",
    "echo \"Current job name  : $JOB_NAME\"\n",
    "echo \"Task index number : $TASK_ID\"\n",
    "echo \"==========================================================\"\n",
    "\n",
    "module load R/4.4.3\n",
    "\n",
    "Rscript ", r_script_path, "\n",
    "\n",
    "End_Time=$(date +\"%s\")\n",
    "echo \"==========================================================\"\n",
    "echo \"Ending on         : $(date)\"\n",
    "Elapsed_Time=$(($End_Time - $Start_Time))\n",
    "echo \"Elapsed time (s)  : $Elapsed_Time\"\n",
    "echo \"==========================================================\"\n"
  )
  writeLines(qsub_content, con = qsub_file_path)

}