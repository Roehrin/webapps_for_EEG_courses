const scriptTag = document.currentScript;
const caller = scriptTag.dataset.caller;
console.log(caller);
const url1 = "docs.html"

function toggleIframeDisplay() {
      const iframeSection = document.getElementById('iframeSection');
      const toggleBtn = document.getElementById('-ToggleIframe');

      console.log('Toggle function called'); // Debug log

      if (iframeSection.style.display === 'none' || iframeSection.style.display === '') {
        iframeSection.style.display = 'block';
		changeIframeUrl(url1.concat("#",caller));
        toggleBtn.textContent = 'Hide Documentation';
        //toggleBtn.classList.add('dropdown-item-checked');
      } else {
        iframeSection.style.display = 'none';
        toggleBtn.textContent = 'Show Documentation';
        //toggleBtn.classList.remove('dropdown-item-checked');
      }
    }

    // Function to change iframe URL (you can call this to load different pages)
    function changeIframeUrl(url) {
      document.getElementById('documentationFrame').src = url;
    }