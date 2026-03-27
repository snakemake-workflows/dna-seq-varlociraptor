function wrapText(value) {
  const max = 80;
  

  let out = '';
  for (let i = 0; i < value.length;) {
    let j = value.indexOf(' ', i + max);
    if (j === -1) {
      out += value.slice(i);
      break;
    }
    out += value.slice(i, j) + '<br>';
    i = j + 1;
  }

  out = out.replace(/\b\d{8}\b/g, id =>
    `<a href='https://pubmed.ncbi.nlm.nih.gov/${id}/' target='_blank'>${id}</a>`
  );

  return `<div style="line-height:1.1 !important">${out}</div>`;
}
