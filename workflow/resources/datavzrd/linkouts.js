function hgvsc_dropdown(value, row) {
  function varsome_link(row) {
    let ref = row["reference allele"];
    let alt = row["alternative allele"];
    if ( alt && ref) {
      let chr = row.chromosome;
      let pos = row.position;
      let build = ( row.build === "GRCh38" ) ? "hg38/" : "hg19/";
      const url = "https://varsome.com/variant/";
      let aa_regex = /^[ACGTacgt]+$/;
      if ( ref.match(aa_regex) && alt.match(aa_regex)) {
        let descriptor = `${chr}:${pos}:${ref}:${alt}`
        return `${url}${build}${descriptor}`;
      }
    }
    return ``
  }

  function genomenexus_link(row) {
      let hgvsg = row.hgvsg
      if ( hgvsg ) {
        let build = ( row.build === "GRCh38" ) ? "grch38." : "";
        const url_suffix = "genomenexus.org/variant/"
        return `https://${build}${url_suffix}${hgvsg}`
      }
      return ``
    }

  function alphamissense_link(row) {
      if ( row.swissprot && row["protein position"] ) {
        let swissprot = row.swissprot.split('.')[0]
        let residue = row["protein position"]
        const url_prefix = "https://alphamissense.hegelab.org/hotspot?uid="
        return `${url_prefix}${swissprot}&resi=${residue}`
      }
      return ``
    }

  let hgvs_regex = /^(.+?):(.+)$/;
  // discard protein id
  let hgvsc = ( value ) ? value.match(hgvs_regex)[2].replace('%3D', '=') : "-";
  let links = "";

  let varsome_url = varsome_link(row);
  if (varsome_url != "") {
      links = `<a class="dropdown-item" href='${varsome_url}' target='_blank' rel='noopener noreferrer' >VarSome</a>`
  }

  let genomenexus_url = genomenexus_link(row)
  if (genomenexus_url != "") {
      links = `${links}<a class="dropdown-item" href='${genomenexus_url}' target='_blank' rel='noopener noreferrer' >Genome Nexus</a>`
  }

  let alphamissense_url = alphamissense_link(row);
  if (alphamissense_url != "") {
      links = `${links}<a class="dropdown-item" href='${alphamissense_url}' target='_blank' rel='noopener noreferrer' >AlphaMissense</a>`
  }
  if ( links ) {
    let dropdown = `<div class="btn-group">
    <button class="btn btn-outline-secondary btn-table btn-sm dropdown-toggle" type="button" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
    ${hgvsc}
    </button>
    <div class="dropdown-menu">
    ${links}
    </div>
    </div>
    `;
    return dropdown;
  }
  return ``
}