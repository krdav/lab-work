# Runtime statistics for Seahawks IncuCyte
Tools for calculating tool runtime/reserved time per Seahawks user. 


### Content
input/
: seahawks_device-log.txt : Log, exported from the IncuCyte software under the "Status" tab and "Logs" sub-tab. Click "Export log" and click "Yes" to exporting the entire log.
: seahawks_vessel-ids.txt : Vessel IDs, exported from the IncuCyte software under the "View" tab. Click "Copy grid to clipboard" and paste text into notepad and save.

output_excel/
: scan_data.xlsx : Data for each scan i.e. each row is one scan with associatied information
: vessel_time.xlsx : Scan data aggregated based on vessel ID i.e. each row is a vessel
: user_time.xlsx : Scan data aggregated based on username (called "Owner") i.e. each row is a user

Column description for user_time.xlsx:
* User reserved time [column: "Reserved time"]. Calculated as the sum of reserved time per vessels. Reversed time for each vessel is calculated by taking the time from start of first scan to end of last scan.
* User scan time [column: "Scan time X fraction"]. Calculated as the sum of scan time per vessels. Scan time for each vessel is calculated by summing all the time spent on all scans for the vessel.


output_plots/
: acc-resrv-user.pdf : Barplot of accumulated user reserved time.
: acc-scan-user.pdf : Barplot of accumulated user scan time.
