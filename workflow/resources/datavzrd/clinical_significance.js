function clin_sig(value) {
  const colorMapping = {
    "unknown": "#6c757d", // gray
    "not_provided": "#343a40", // dark gray
    "other": "#adb5bd", // light gray
    "benign": "#218838", // green
    "benign/likely_benign": "#218838", // dark green
    "likely_benign": "#28a745", // green
    "protective": "#20c997", // teal
    "uncertain_significance": "#ffc107", // yellow
    "conflicting_interpretations_of_pathogenicity": "#17a2b8", // cyan
    "association": "#007bff", // blue
    "affects": "#6f42c1", // purple
    "drug_response": "#6610f2", // dark purple
    "risk_factor": "#e3342f", // red
    "likely_pathogenic": "#dc3545", // dark red
    "pathogenic/likely_pathogenic": "#bd2130", // deep red
    "pathogenic": "#bd2130" // dark red
  };
  
  const splitValues = value.split(",").map(item => `<span style="color: white; background-color: ${colorMapping[item.trim()]}"class="badge">${item.trim()}</span>`);
  return splitValues.join(' ');  
}
